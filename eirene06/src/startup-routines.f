C EIRENE06 COMPILATION
C ===== SOURCE: areaa.f
C  OCT 03: return best guess, even if R out of range, rather then exit
C
      FUNCTION AREAA (R,N,ARCA,Y,EP1R,ELLR)
C
C   FOR LEVGEO=2 OPTION
C
C   THIS FUNCTION INTERPOLATES LINEAR STANDARD MESH PARAMETERS ON POINTS
C   LYING BETWEEN THE SURFACES OF THIS MESH
C   AND EVALUATES THE ENCLOSED (CROSS SECTIONAL) AREA AND THE ARCLENGTH
C
C  INPUT:
C     R
C     N MUST BE GIVEN SUCH THAT RSURF(N)<=R<=RSURF(N+1)
C  OUTPUT:
C     AREAA(R)=AREA WITHIN EIRENE-SURFACE LABELED BY R (RSURF)
C                   I.E. CROSS SECTIONAL AREA (CM**2)
C     ARCA(R)=ARCLENGTH AT R
C     Y(R)
C     EP1R(R)
C     ELLR(R)
C     TRIA(R) (TO BE WRITTEN)
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: R
      REAL(DP), INTENT(OUT) :: ARCA, Y, EP1R, ELLR
      INTEGER, INTENT(IN) :: N
      REAL(DP) :: RRN, APB, RRI, RRD, XL1QQ, Q1, XL1, XL1Q, AREAA
      INTEGER :: NP
C
      NP=N+1
      IF (N.LE.0.OR.N.GE.NR1STM) GOTO 998
      IF (R.LT.RSURF(N).OR.R.GT.RSURF(NP)) GOTO 999
      RRI=RSURF(N)
      RRD=RSURF(NP)-RRI
      RRN=(R-RRI)/RRD
C
      ELLR=ELL(N)+RRN*(ELL(NP)-ELL(N))
      EP1R=EP1(N)+RRN*(EP1(NP)-EP1(N))
 100  Y=ELLR*R
C
      APB=R+Y+EPS60
      XL1=(R-Y)/APB
      XL1Q=XL1*XL1
      XL1QQ=XL1Q*XL1Q
      Q1=(16.-0.75*XL1QQ)/(64.-16.*XL1Q)
C
      AREAA=PIA*Y*R
      ARCA=APB*Q1*PIA*4.
C
      RETURN
C
998   WRITE (iunout,*) 'ERROR IN FUNCTION AREAA: N, NR1STM ',N,NR1STM
      CALL EXIT_OWN(1)
999   WRITE (iunout,*) 'ERROR IN FUNCTION AREAA: R,N ',R,N
      IF (R.LT.RSURF(N)) THEN
        ELLR=ELL(N)
        EP1R=EP1(N)
      ELSE
        ELLR=ELL(NP)
        EP1R=EP1(NP)
      ENDIF
      GOTO 100
      END
C ===== SOURCE: ergod.f
C
      SUBROUTINE ERGOD

C  MODIFICATIONS TO INPUT, IN ORDER TO PERFORM A CELL VOLUME ESTIMATION RUN
C  BASED UPON AN ERGODIC PRINCIPLE

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CTRCEI
      USE COMPRT
      USE COMNNL
      USE COMSOU
      USE CTEXT
      USE CSDVI
      USE CLGIN
      USE COMXS

      IMPLICIT NONE

      REAL(DP) :: SORL4
      INTEGER :: ISTR, ISPR, ISTR2, J, ISTR_ERG
      INTEGER, EXTERNAL :: IDEZ
C
C  1.)  MAKE ALL NON-TRANSPARENT SURFACES 100% REFLECTING, 1.5 EV MONOENERGETIC
C       "THERMAL RE-EMISSION" MODEL
C
      DO J=1,NLIMPS
        IF (ILIIN(J).GT.0) THEN
          ILIIN(J)=1
          ILSPT(J)=0
          ILSIDE(J)=0
          EWALL(J)=1.5
C         RINTEG(J)=0.0
          DO ISPZ=1,NSPZ
            EXPIL(ISPZ,J)=0.0
C  THESE SETTINGS: EITHER RESULT ON PDENA(1,...), THEN: PDENM(1,....)=0.
C                  OR     RESULT ON PDENM(1,...), THEN: PDENA(1,....)=0.
            TRANSP(ISPZ,1,J)=0.0
            TRANSP(ISPZ,2,J)=0.0
            RECYCF(ISPZ,J)=0.0
            RECYCT(ISPZ,J)=1.0
            ISRF(ISPZ,J)=0
            ISRT(ISPZ,J)=1
          ENDDO
        ENDIF
      ENDDO
C
C  2.)  TURN OFF ALL VOLUME PROCESSES
C
      DO IATM=1,NATM
        NRCA(IATM)=-1
      ENDDO
      DO IMOL=1,NMOL
        NRCM(IMOL)=-1
      ENDDO
      DO IION=1,NION
        NRCI(IION)=-1
      ENDDO
C
C  3.)  TURN ON TIME-HORIZON, CENSUS ARRAYS, ETC.
c       CONSIDER ONLY ONE STRATUM: ISTR_ERG
C       PLUS ONE FURTHER STRATUM FOR CENSUS
C
      DO ISTR=1,NSTRAI-1
        IF (FLUX(ISTR).GT.0.AND.NPTS(ISTR).GT.0) THEN
          ISTR_ERG=ISTR
C   SET 4TH DIGIT OF SORLIM = 1
          DO J=1,NSRFSI(ISTR_ERG)
            SORL4=IDEZ(INT(SORLIM(J,ISTR_ERG)),4,4)
            SORLIM(J,ISTR_ERG)=SORLIM(J,ISTR_ERG)-SORL4*1000+1000
          ENDDO
C  TURN OFF ALL FURTHER STRATA
          DO ISTR2=ISTR+1,NSTRAI-1
            NPTS(ISTR2)=0
            FLUX(ISTR2)=0.
          ENDDO
        ENDIF
      ENDDO
      NPTS(NSTRAI)=0
      FLUX(NSTRAI)=0
C
C  4.)  TURN ON PRINTOUT OF PARTICLE DENSITY, VOLUME AND CENSUS FLUXES
C
C     TRCHST=.FALSE.
      TRCAMD=.FALSE.
      DO ISPR=1,NSURPR
        IF (NPRSRF(ISPR).EQ.NLIM+NSTSI) GOTO 401
      ENDDO
      NSURPR=NSURPR+1
      NPRSRF(NSURPR)=NLIM+NSTSI
401   CONTINUE
C
      NVOLPR=4
C  PRINT CELL VOLUMES
      NPRTLV(1)=-NTALO
      NFLAGV(1)=3
      NSPEZV(1,1)=1
      NSPEZV(1,2)=1
      NTLVFL(1)=91
C  PRINT EITHER THERMAL ATOMS, SPECIES 1
      NPRTLV(2)=1
      NFLAGV(2)=3
      NSPEZV(2,1)=1
      NSPEZV(2,2)=1
      NTLVFL(2)=92
C  OR PRINT THERMAL MOLECULES, SPECIES 1
      NPRTLV(3)=2
      NFLAGV(3)=3
      NSPEZV(3,1)=1
      NSPEZV(3,2)=1
      NTLVFL(3)=93
C  OR PRINT USER DEFINED OUTPUT FROM TALUSR
      NPRTLV(4)=0
      NFLAGV(4)=3
      NSPEZV(4,1)=1
      NSPEZV(4,2)=1
      NTLVFL(4)=99
C
C  TENTATIVELY: TURN ON STATISTICAL ERROR ESTIMATOR FOR ATOMS, SPECIES 1
C               CORRECT THAT IN CASE OF MOLECULES, SPECIES 1, LATER
      NSIGVI=1
      NSIGCI=0
      NSIGI =NSIGVI+NSIGSI+NSIGCI
      IIH(1)=1
C  ...OR: IIH(1)=2, SEE: PARTICLE LOOP IN SUBR. MCARLO, AFTER FIRST HISTORY
      IGH(1)=1
C
      RETURN
      END
C ===== SOURCE: find_param.f
C
!pb  11.12.06:  allow letters 'f' or 't' in case name of fem or tetrahedron 
!pb             calculation
C
      SUBROUTINE FIND_PARAM
C
C   SEARCH FOR PARAMETERS AND SET DEFAULT VALUES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER :: INDGRD(3), INDPRO(12)
      INTEGER, ALLOCATABLE :: INDSRC(:)
      INTEGER :: ISTRA, NSTRAI, NFLR, ISOR, NSRFSI, NRADD,
     .           NREACI, NSTSI, NLIMI, NVOLPL, NSP, ICO,
     .           NSURPR, NVOLPR, NPRNLI, NCHORI,
     .           NCHENI, NSIGI_BGK, NSIGSI, ID, NSIGVI,
     .           NSIGI_COP, NR1ST, NRSEP, NTIME0,
     .           NP1, NP2, NRKNOT, NRPLG, NPPLG, IERROR, IUNIN,
     .           NITER0, K, NTPER, NTTRA, NCOOR, NTET,
     .           NT3RD, NTSEP, NTRII, NP2ND, I, J, NPPER, NPSEP, NPPLA,
     .           NSIGCI, IREAD, NCOPI, NCOPII, NCOPIE, NREAC_ADD,
     .           IPLS, NRC, NRE, NLINES, I2, LL, NB1, NB2, NB3, NS1, 
     .           NS2, NS3, INM1, INM2, INM3, INMDL
      REAL(DP) :: SORIND, SORLIM, DUMM1, ROA, ZAA, ZZA, ZGA, YAA, YYA,
     .            ZIA, YP, XP, YIA, YGA
      LOGICAL :: NLSCL, NLTEST, NLANA, NLDRFT, NLCRR, NLERG, LIDENT,
     .           LHABER
      LOGICAL :: NLSLB, NLCRC,  NLELL, NLTRI,  NLPLG, NLFEM, NLTET,
     .           NLGEN
      LOGICAL :: NLRAD, NLPOL,  NLTOR, NLADD,  NLMLT, NLTRIM
      LOGICAL :: NLTRA, NLTRT, NLTRZ
      LOGICAL :: PLTL2D, PLTL3D, LRPSCUT
      CHARACTER(72) :: ZEILE, CASENAME, FILENAME, ULINE
C
C  SET DEFAULT VALUES FOR PARAMETERS
C
C  GEOMETRY
      N1ST=1
      N2ND=1
      N3RD=1
!pb      NADD=1
      NADD=0
      NTOR=1
      NRTAL=0
      NLIM=1
      NSTS=1
      NPLG=1
      NPPART=1
      NKNOT=1
      NTRI=1
      NTRII=0
      NTETRA=1
      NCOORD=1
      NOPTIM=1
      NOPTM1=1
C  PRIMARY SOURCE
      NSTRA=1
      NSRFS=1
      NSTEP=1
      NGITT=0
C  SPECIES AND TALLIES
      NATM=0
      NMOL=0
      NION=0
      NPLS=1
      NPHOT=0
      NADV=0
      NADS=0
      NCLV=0
      NSNV=0
      NALV=0
      NALS=0
      NAIN=1
      NCOP=0
      nbgk=0
      NADSPC=0
      NPLSTI=1
      NPLSV=1
C  STATISTICS
      NSD=1
      NSDW=1
      NCV=1
C  ATOMIC DATA
      NREAC=1
      NREAC_ADD=0
      NREC=1
      NRDS=1
      NRCX=1
      NREL=1
      NRPI=1
      NPTRGT=1

      NCHOR=0
      NCHEN=0
      NPRNL=0

C  Stellarator geometry?
C  NGEOM_USR = 1  ==> Stellarator-Geometrie
C  NGEOM_USR = 0  ==> keine Stellarator-Geometrie
      NGEOM_USR=0

C  Input from coupling routine
C  NCOUP_INPUT = 0  ==> NO PLASMA INPUT FROM COUPLING ROUTINE
C  NCOUP_INPUT = 1  ==> PLASMA INPUT FROM COUPLING ROUTINE
      NCOUP_INPUT = 1

C  Sum over Strata
C  NSMSTRA = 0  ==> SUM OVER STRATA IS NOT PERFORMED
C  NSMSTRA = 1  ==> SUM OVER STRATA IS PERFORMED
      NSMSTRA=1
C
C  CALULATION OR STORAGE OF ATOMIC DATA
C  NSTORAM=0,9.    =0: MINIMUM STORAGE, MAXIMUM CALULATION
C                  =9: MAXIMUM STORAGE, MINIMUM CALULATION
      NSTORAM=9
C
C  SPATIALLY RESOLVED SURFACE TALLIES?
C  NGSTAL = 0  ==> NO SPACE FOR SPATIALLY RESOLVED SURFACE TALLIES
C  NGSTAL = 1  ==> SPATIALLY RESOLVED SURFACE TALLIES ARE COMPUTED
      NGSTAL=0

C  NUMBER OF BACKGROUND SPECTRA
      NBACK_SPEC=0
C
C  INITIALIZE SOME DATA AND SET DEFAULTS
C

C
C  UNIT NUMBER FOR INPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15
      IUNIN=1

      REWIND IUNIN
C
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1).EQ.'*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NMACH,NMODE,NTCPU,NFILE,NITER0,NITER,
     .                    NTIME0,NTIME

      READ (IUNIN,'(A72)') ZEILE
      IF ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') + INDEX(ZEILE,'T') +
     .     INDEX(ZEILE,'t')) == 0) THEN
        READ (ZEILE,6666) NOPTIM,NOPTM1,NGEOM_USR,NCOUP_INPUT,
     .                    NSMSTRA,NSTORAM,NGSTAL,NRTAL,NREAC_ADD
        READ (IUNIN,6665) NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
     .                    NLERG,LIDENT,LHABER
      ELSE
        READ (ZEILE,6665) NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
     .                    NLERG,LIDENT,LHABER
      END IF

      NSTORAM = MIN(NSTORAM,9)
      IF (NSTORAM < 9) NSTORAM = 0

      WRITE (iunout,*) '*** 1. DATA FOR OPERATING MODE '

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C
C  READ DATA FOR STANDARD MESH, 200---299
C
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) (INDGRD(J),J=1,3)
C
C INPUT SUB-BLOCK 2A
C
C  RADIAL MESH
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLRAD
      IREAD=0
      IF (NLRAD) THEN
C
        READ (IUNIN,6665) NLSLB,NLCRC,NLELL,NLTRI,NLPLG,
     .                    NLFEM,NLTET,NLGEN
c slmod begin
        IF (NLTET) THEN 
          READ (IUNIN,6667) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
        ELSE
          READ (IUNIN,6666) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
        ENDIF
c
c        READ (IUNIN,6666) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
c slmod end
        N1ST = MAX(N1ST,NR1ST)
        IF (NLPLG) NPLG = MAX(NPLG,NRPLG)
        IF (NLPLG) NPPART = MAX(NPPART,NPPLG)
        NKNOT = MAX(NKNOT,NRKNOT)
        NCOORD = MAX(NCOORD,NCOOR)
        IF (INDGRD(1).LE.5) THEN
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            READ (IUNIN,*)
            IF (NLELL.OR.NLTRI) THEN
              READ (IUNIN,*)
              READ (IUNIN,*)
              IF (NLTRI) THEN
                READ (IUNIN,*)
              ENDIF
            ENDIF
          ENDIF
          IF (NLPLG) THEN
            READ (IUNIN,*)
            READ (IUNIN,6666) (NP1,NP2,K=1,NPPLG)
            DO 212 I=1,NR1ST
              READ (IUNIN,6664) (XP,YP,    J=1,NRPLG)
212         CONTINUE
          ENDIF
          IF (NLFEM) THEN
            NTRII=NR1ST
            NTRI=MAX(NTRI,NR1ST)
            READ (IUNIN,'(A72)') ZEILE
            IF ((INDEX(ZEILE,'CASE')==0) .AND. 
     .          (INDEX(ZEILE,'case') == 0)) THEN
              DO WHILE ((ZEILE(1:1) .NE. '*') .AND.
     .                  ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') +
     .                    INDEX(ZEILE,'T') + INDEX(ZEILE,'t')) == 0))
                READ (IUNIN,'(A72)') ZEILE
              END DO
              IREAD=1
            ELSE
              READ (ZEILE(6:),'(A66)') CASENAME
              CASENAME=ADJUSTL(CASENAME)
              I2=INDEX(CASENAME,' ')
              LL=LEN_TRIM(CASENAME)

              FILENAME=CASENAME(1:LL) // '.npco_char'
              OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .              FORM='FORMATTED')
              ZEILE='*   '      
              DO WHILE (ZEILE(1:1) == '*')
                 READ (30,'(A)') ZEILE
              END DO

              READ (ZEILE,*) NRKNOT
              CLOSE (UNIT=30)
              NKNOT=NRKNOT

              FILENAME=CASENAME(1:LL) // '.elemente'
              OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .              FORM='FORMATTED')

              ZEILE='*   '      
              DO WHILE (ZEILE(1:1) == '*')
                 READ (30,'(A)') ZEILE
              END DO

              READ (ZEILE,*) NTRII
              CLOSE (UNIT=30)
              NTRI=NTRII+1
              NR1ST=NTRI

              FILENAME=CASENAME(1:LL) // '.neighbors'
              OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .              FORM='FORMATTED')
              
              ZEILE='*   '
              DO WHILE (ZEILE(1:1) == '*')
                 READ (30,'(A100)') ZEILE
              END DO
              
              NGITT=1
              DO I=1,NTRII
                 READ (30,*) ID, NB1, NS1, INM1,
     .                           NB2, NS2, INM2, 
     .                           NB3, NS3, INM3
                 IF (INM1 /= 0) NGITT = NGITT + 1
                 IF (INM2 /= 0) NGITT = NGITT + 1
                 IF (INM3 /= 0) NGITT = NGITT + 1
              END DO
              
              CLOSE (UNIT=30)
              
            END IF
          ENDIF
          IF (NLTET) THEN
            NTET=NR1ST
            NTETRA=MAX(NTETRA,NR1ST)
            READ (IUNIN,'(A72)') ZEILE
            ULINE = ZEILE
            CALL UPPERCASE(ULINE)
            DO WHILE ((INDEX(ULINE,'CASE') > 0) .OR. 
     .                ((ULINE(1:1) .NE. '*') .AND.
     .                 ((INDEX(ULINE,'F') + INDEX(ULINE,'T')) == 0)))
              READ (IUNIN,'(A72)') ZEILE
              ULINE = ZEILE
              CALL UPPERCASE(ULINE)
            END DO
            IREAD=1
          ENDIF
        ELSEIF (INDGRD(1).EQ.6) THEN
C  IS THERE ONE MORE LINE, OR IS NLPOL THE NEXT VARIABLE
          READ (IUNIN,'(A72)') ZEILE
          DO WHILE ((ZEILE(1:1) .NE. '*') .AND.
     .         ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') +
     .         INDEX(ZEILE,'T') + INDEX(ZEILE,'t')) == 0))
             READ (IUNIN,'(A72)') ZEILE
          END DO
          IREAD=1
        ENDIF
      ENDIF
      NTRIS=NTRII
      NKNOTS=NKNOT
C
C  POLOIDAL MESH
C
C INPUT SUB-BLOCK 2B
C
      IF (IREAD == 0) READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLPOL
C
      READ (IUNIN,*)
      READ (IUNIN,6666) NP2ND,NPSEP,NPPLA,NPPER
      IF (INDGRD(2).LE.5) THEN
        READ (IUNIN,6664) YIA,YGA,YAA,YYA
      ENDIF
      N2ND = MAX(N2ND,NP2ND)
C
C  TOROIDAL MESH
C
C INPUT SUB-BLOCK 2C
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLTOR
C
      READ (IUNIN,6665) NLTRZ,NLTRA,NLTRT
      READ (IUNIN,6666) NT3RD,NTSEP,NTTRA,NTPER
      IF (INDGRD(3).LE.5) THEN
        READ (IUNIN,6664) ZIA,ZGA,ZAA,ZZA,ROA
      ELSEIF (INDGRD(3).EQ.6) THEN
      ENDIF
      IF (NLTOR.AND.NLTRA) NTTRA=NT3RD
      N3RD = MAX(N3RD,NT3RD)
      NTOR = MAX(NTOR,NTTRA)
C
C  MESH MULTIPLICATION
C
C INPUT SUB-BLOCK 2D
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLMLT
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE ((ZEILE(1:1) .NE. '*') .AND.
     .          ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') +
     .            INDEX(ZEILE,'T') + INDEX(ZEILE,'t')) == 0))
         READ (IUNIN,'(A72)') ZEILE
      END DO
      IF (ZEILE(1:1) .EQ. '*') READ (IUNIN,'(A72)') ZEILE
C
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
C
      READ (ZEILE,6665) NLADD
C
      IF (NLADD) THEN
        READ (IUNIN,6666) NRADD
        NADD = MAX(NADD,NRADD)
      ENDIF

C  FIND START OF NEXT INPUT BLOCK: 3A

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C
C  READING FOR INPUT BLOCK 2 DONE
C
C
C  READ DATA FOR NON DEFAULT SURFACE MODELS ON STANDARD SURFACES
C  300--349
C
      WRITE (iunout,*) '*** 3A. DATA FOR NON DEFAULT STANDARD SURFACES'
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ(ZEILE,6666) NSTSI
      IF (NTIME.GE.1) NSTSI = NSTSI + 1
      NSTS = MAX(NSTS,NSTSI)

C  FIND START OF NEXT INPUT BLOCK: 3B

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C     READ DATA FOR ADDITIONAL SURFACES 350--399
C
      WRITE (iunout,*) '*** 3B. DATA FOR ADDITIONAL SURFACES           '
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NLIMI
      NLIM = MAX(NLIM,NLIMI)

C  FIND START OF NEXT INPUT BLOCK: 4

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR SPECIES SPECIFICATION AND ATOMIC PHYSICS MODULE
C  400--499
C
      WRITE (iunout,*) '*** 4. DATA FOR SPECIES SPECIFICATION AND   '
      WRITE (iunout,*) '       ATOMIC PHYSICS MODULE                '
C
      READ (IUNIN,*)
      WRITE (iunout,*) 
     .  '       ATOMIC REACTION CARDS, NREACI DATA FIELDS'
      READ (IUNIN,*) NREACI
      NREAC = MAX(NREAC,NREACI)+NREAC_ADD
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
      WRITE (iunout,*) 
     .  '*4A.   NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES'
      READ (IUNIN,*) NATMI
      NATM = MAX(NATM,NATMI)

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ NEUTRAL MOLECULES SPECIES CARDS
C
      WRITE (iunout,*) 
     .  '*4B.   NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES'
      READ (IUNIN,*) NMOLI
      NMOL = MAX(NMOL,NMOLI)

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ TEST PARTICLE IONS SPECIES CARDS
C
      WRITE (iunout,*) '*4C.   TEST IONS SPECIES CARDS, NIONI SPECIES '
      READ (IUNIN,*) NIONI
      NION = MAX(NION,NIONI)

C  FIND START OF NEXT INPUT BLOCK: 4

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO

C  FIND START OF NEXT INPUT BLOCK: 4D
      NPHOTI=0
      IF (ZEILE(1:3) == '***') GOTO 500
      WRITE (iunout,*) '*4D.   PHOTONS SPECIES CARDS, NPHOTI SPECIES '
      READ (IUNIN,*) NPHOTI
      NPHOT = MAX(NPHOT,NPHOTI)

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR PLASMA-BACKGROUND , 500--599
C
      WRITE (iunout,*) '*** 5. DATA FOR PLASMA BACKGROUND            '
C
C  READ BULK IONS SPECIES CARDS
C
      WRITE (iunout,*) '*5A.   BULK ION SPECIES CARDS, NPLSI SPECIES '

      READ (IUNIN,'(A72)') ZEILE
 500  DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ(ZEILE,6666) NPLSI
      NPLS = MAX(NPLS,NPLSI)

      ICO = 0
      DO IPLS=1,NPLSI
        READ (IUNIN,'(A72)') ZEILE
        ULINE = ZEILE
        CALL UPPERCASE(ULINE)
        INMDL=INDEX(ULINE,'FORT')+INDEX(ULINE,'SAHA')+
     .        INDEX(ULINE,'CORONA')+
     .        INDEX(ULINE,'BOLTZMANN')+INDEX(ULINE,'COLRAD')
        IF (INMDL > 0) ICO = ICO + 1
        READ (ZEILE(33:35),'(I3)') NRC
        DO K=1,NRC
          READ (IUNIN,*)
          READ (IUNIN,*)
        END DO
!pb        IF (VERIFY(ZEILE(57:66),' ') > 0) THEN
        IF (INMDL > 0) THEN
          NRE=0
          IF (VERIFY(ULINE(INMDL+11:),' ') > 0) 
     .      READ(ZEILE(INMDL+11:),*) NRE
          NRE=MAX(NRE,1)
          NLINES=1
          IF (INDEX(ULINE(INMDL:),'COLRAD') > 0) NLINES=NRE
          DO K=1,NLINES
             READ (IUNIN,*)
          END DO
        END IF
      END DO
      IF (ICO > 0) NREAC=NREAC+1

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) == '*')
         READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) (INDPRO(J),J=1,12)

      NPLSTI = 1
      IF ((INDPRO(2) < 0) .OR. (INDPRO(2) > 9)) NPLSTI=NPLS
      NPLSV = NPLS
      IF (ABS(INDPRO(4)) > 9) NPLSV = 1

C  FIND START OF NEXT INPUT BLOCK: 6

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ  DATA FOR REFLECTION MODEL  600--699
C
      WRITE (iunout,*) '*** 6. GENERAL DATA FOR REFLECTION MODEL    '
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLTRIM
      NFLR=0
      IF (NLTRIM) THEN
        READ (IUNIN,'(A72)') ZEILE
        IF (INDEX(ZEILE,'PATH')+INDEX(ZEILE,'path').NE.0) THEN
C  PATH SPECIFICATION FOR DATA BASE FOUND
          READ (IUNIN,'(A72)') ZEILE
          DO WHILE ((INDEX(ZEILE,'ON')+INDEX(ZEILE,'on')) > 0)
            NFLR=NFLR+1
            READ (IUNIN,'(A72)') ZEILE
          END DO
        ENDIF
      ENDIF

      IF (NFLR > 0) THEN
        NHD6 = NFLR
      ELSE
        NHD6 = 8
      END IF

C  FIND START OF NEXT INPUT BLOCK: 7

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR PRIMARY SOURCE  700--799
C
      WRITE (iunout,*) '*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA'
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NSTRAI

      NSTEP = 1
      ALLOCATE (INDSRC(NSTRAI))
      READ (IUNIN,6666) (INDSRC(J),J=1,NSTRAI)
      IF (ANY(INDSRC == 6)) THEN
        DO J=NSTRAI,1,-1
          IF (INDSRC(J) == 6) THEN
            NSTEP = J
            EXIT
          END IF
        END DO
      END IF
      READ (IUNIN,*)
      DO ISTRA=1,NSTRAI
        IF (INDSRC(ISTRA) == 6) CYCLE
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,6666) NSRFSI
        NSRFS = MAX(NSRFS,NSRFSI)

        DO I=1,NSRFSI
          READ (IUNIN,*)
          READ (IUNIN,6664) DUMM1, SORLIM, SORIND
          ISOR = INT(SORLIM)
          DO WHILE (ISOR > 0)
            ID = MOD(ISOR,10)
            IF (ID == 4) NSTEP = MAX(NSTEP,INT(SORIND))
            ISOR = ISOR / 10
          END DO
          READ (IUNIN,*)
          READ (IUNIN,*)
        END DO
        READ (IUNIN,*)
        READ (IUNIN,*)
C
      END DO
      DEALLOCATE (INDSRC)
      IF (NTIME.GE.1) NSTRAI = NSTRAI + 1
      NSTRA = MAX(NSTRA,NSTRAI)

      READ (IUNIN,'(A72)') ZEILE
C
C     READ ADDITIONAL DATA FOR SOME SPECIFIC ZONES
C
      WRITE (iunout,*) '*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES '

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO

C
C  READ DATA FOR STATISTICS AND NONANALOG MODEL, 900--999
C
      WRITE (iunout,*) '*** 9. DATA FOR STATISTIC AND NONANALOG MODEL '

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  DATA FOR STANDARD DEVIATION
      WRITE (iunout,*) '       CARDS FOR STANDARD DEVIATION '
      READ (IUNIN,6666) NSIGVI,NSIGSI,NSIGCI,NSIGI_BGK,NSIGI_COP
      NSD = MAX(NSD,NSIGVI)
      NSDW = MAX(NSDW,NSIGSI)
      NCV = MAX(NCV,NSIGCI)

C  FIND START OF NEXT INPUT BLOCK: 7

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C   READ DATA FOR ADDITIONAL AND SURFACE AVERAGED TALLIES
      WRITE (iunout,*) 
     .  '*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION   '
      WRITE (iunout,*) '        ESTIMATORS  AND ALGEBRAIC EXPRESSIONS  '
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NADVI,NCLVI,NALVI,NADSI,NALSI,NADSPC
      NADV = MAX(NADV,NADVI)
      NCLV = MAX(NCLV,NCLVI)
      NALV = MAX(NALV,NALVI)
      NADS = MAX(NADS,NADSI)
      NALS = MAX(NALS,NALSI)
C
      WRITE (iunout,*) '*** 10A. DATA FOR ADDITIONAL TALLIES           '
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10B. DATA FOR COLLISION ESTIMATORS         '
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS     '
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES   '
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES '
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10F. DATA FOR SPECTRA '
      READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3) .NE. '***') THEN
        READ (IUNIN,'(A72)') ZEILE
        DO WHILE (ZEILE(1:3) .NE. '***')
          READ (IUNIN,'(A72)') ZEILE
        END DO
      END IF
C
C   READ DATA FOR NUMERICAL AND GRAPHICAL OUTPUT 1100--1199
C
      WRITE (iunout,*) 
     .  '*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT  '
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
      READ (IUNIN,*)
      READ (IUNIN,6666) NVOLPR
      NVOLPR=MIN(NVOLPR,100)
      DO J=1,NVOLPR
        READ (IUNIN,*)
      END DO
C
      READ (IUNIN,6666) NSURPR
      NSURPR=MIN(NSURPR,100)
      DO J=1,NSURPR
        READ (IUNIN,*)
      END DO
C
      READ (IUNIN,'(A72)') ZEILE
      CALL UPPERCASE(ZEILE)
      DO WHILE (SCAN(ZEILE,'*FT') == 0)
        READ (IUNIN,'(A72)') ZEILE
        CALL UPPERCASE(ZEILE)
      END DO
C
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      READ (ZEILE(16:16),'(L1)') LRPSCUT

      DO J=1,13
        READ (IUNIN,*)
      END DO
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  DATA FOR PLOTS OF VOLUME AVERAGED TALLIES
      READ (ZEILE,6666) NVOLPL
C
C
      NPLT = 1
      IF (NVOLPL > 0) THEN
        READ (IUNIN,*)
        IF (LRPSCUT) READ (IUNIN,*)
        DO J=1,NVOLPL
          READ (IUNIN,'(A72)') ZEILE
          DO WHILE (ZEILE(1:1) .EQ. '*')
            READ (IUNIN,'(A72)') ZEILE
          END DO
          READ (ZEILE,6666) NSP
          NPLT = MAX(NPLT, NSP)
          READ (IUNIN,6665) PLTL2D,PLTL3D
          READ (IUNIN,*)
          IF (PLTL2D) THEN
            READ (IUNIN,*)
            DO I=1,NSP
              READ (IUNIN,*)
            END DO
          ENDIF
          IF (PLTL3D) THEN
            READ (IUNIN,*)
            READ (IUNIN,*)
            DO I=1,NSP
              READ (IUNIN,*)
            END DO
            READ (IUNIN,*)
          ENDIF
        END DO
      END IF
C
C  SKIP INPUT LINES, UNTIL INPUT BLOCK 12 STARTS
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR DIAGNOSTIC MODULE  1200--1299
C
      WRITE (iunout,*) '*** 12. DATA FOR DIAGNOSTIC MODULE '
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NCHORI,NCHENI
      NCHOR = MAX(NCHOR,NCHORI)
      NCHEN = MAX(NCHEN,NCHENI)
      IF (NCHORI > 0) THEN
        NREAC=NREAC+1
        NADV=NADV+10
      END IF

C  SKIP READING REST OF THIS BLOCK
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR NONLINEAR MODE  1300--1399
C
      WRITE (iunout,*) 
     .  '*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION '
C
      READ (IUNIN,6666) NPRNLI
      IF (NLERG.AND.NPRNLI.LE.0) THEN
C  NO TIME HORIZON DEFINED, DESPITE NLERG=.TRUE.
C  THEREFORE: SET A DEFAULT TIME HORIZON HERE
        NPRNLI=100
      ENDIF
      NPRNL = MAX(NPRNL,NPRNLI)

      IF (NPRNLI > 0) THEN
        READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1).NE.'*') THEN
          READ (IUNIN,*)
C   READ DATA FOR SNAPSHOT TALLIES
          READ (IUNIN,*)
          WRITE (iunout,*) '*** 13A. DATA FOR SNAPSHOT TALLIES'
          READ (IUNIN,6666) NSNVI
          NSNV = MAX(NSNV,NSNVI)
        END IF
      END IF
C  SKIP READING REST OF THIS BLOCK
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  INPUT BLOCK 14 BEGIN
C
C  READ DATA IN INTERFACING SUBROUTINE INFCOP  1400 -- 1499
C
      WRITE (iunout,*) '*** 14. DATA FOR INTERFACING ROUTINE "INFCOP" '
      IF (NMODE.EQ.0) THEN
        READ (IUNIN,6666) NAINI,NCOPII,NCOPIE
        NCOPI=NCOPIE
        NAIN = MAX(NAIN,NAINI)
        NCOP = MAX(NCOP,NCOPI)
        NCPV = NCOP
      ELSE
        CALL IF0PRM(IUNIN)
      ENDIF

      REWIND IUNIN
C
      WRITE (iunout,*) 'N1ST = ',N1ST
      WRITE (iunout,*) 'N2ND = ',N2ND
      WRITE (iunout,*) 'N3RD = ',N3RD
      WRITE (iunout,*) 'NADD = ',NADD
      WRITE (iunout,*) 'NTOR = ',NTOR
      WRITE (iunout,*) 'NRTAL = ',NRTAL
      WRITE (iunout,*) 'NLIM = ',NLIM
      WRITE (iunout,*) 'NSTS = ',NSTS
      WRITE (iunout,*) 'NPLG = ',NPLG
      WRITE (iunout,*) 'NPPART = ',NPPART
      WRITE (iunout,*) 'NKNOT = ',NKNOT
      WRITE (iunout,*) 'NTRI = ',NTRI
      WRITE (iunout,*) 'NTETRA = ',NTETRA
      WRITE (iunout,*) 'NCOORD = ',NCOORD
      WRITE (iunout,*) 'NOPTIM = ',NOPTIM
      WRITE (iunout,*) 'NOPTM1 = ',NOPTM1
      WRITE (iunout,*) 'NSTRA = ',NSTRA
      WRITE (iunout,*) 'NSRFS = ',NSRFS
      WRITE (iunout,*) 'NSTEP = ',NSTEP
      WRITE (iunout,*) 'NATM = ',NATM
      WRITE (iunout,*) 'NMOL = ',NMOL
      WRITE (iunout,*) 'NION = ',NION
      WRITE (iunout,*) 'NPLS = ',NPLS
      WRITE (iunout,*) 'NADV = ',NADV
      WRITE (iunout,*) 'NADS = ',NADS
      WRITE (iunout,*) 'NCLV = ',NCLV
      WRITE (iunout,*) 'NSNV = ',NSNV
      WRITE (iunout,*) 'NALV = ',NALV
      WRITE (iunout,*) 'NALS = ',NALS
      WRITE (iunout,*) 'NAIN = ',NAIN
      WRITE (iunout,*) 'NCOP = ',NCOP
      WRITE (iunout,*) 'NBGK = ',NBGK
      WRITE (iunout,*) 'NSD = ',NSD
      WRITE (iunout,*) 'NSDW = ',NSDW
      WRITE (iunout,*) 'NCV = ',NCV
      WRITE (iunout,*) 'NREAC = ',NREAC
      WRITE (iunout,*) 'NREC = ',NREC
      WRITE (iunout,*) 'NRDS = ',NRDS
      WRITE (iunout,*) 'NRCX = ',NRCX
      WRITE (iunout,*) 'NREL = ',NREL
      WRITE (iunout,*) 'NRPI = ',NRPI
      WRITE (iunout,*) 'NGEOM_USR = ',NGEOM_USR
      WRITE (iunout,*) 'NCOUP_INPUT = ',NCOUP_INPUT
      WRITE (iunout,*) 'NSMSTRA = ',NSMSTRA
      WRITE (iunout,*) 'NSTORAM = ',NSTORAM
      WRITE (iunout,*) 'NGSTAL = ',NGSTAL
      WRITE (iunout,*) 'NRPES  = ',NRPES
C
      RETURN
C
6664  FORMAT (6E12.4)
6665  FORMAT (12(5L1,1X))
6666  FORMAT (12I6)
c slmod begin
6667  FORMAT (12I8)
c slmod end

      END
C ===== SOURCE: grid_1.f
C
C
      SUBROUTINE GRID_1(C,NT,N2,N3,C1,C2,C3,CT,IL)
C
C  MAKE A 1D GRID C(I),I=1,NT
C  CONSISTS OF UP TO 3 SECTIONS OF EQUIDISTANT GRID POINTS C
C  C(1) =C1
C  C(N2)=C2
C  C(N3)=C3
C  C(NT)=CT
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: C1, C2, C3, CT
      REAL(DP), INTENT(OUT) :: C(*)
      INTEGER, INTENT(IN) :: NT, IL
      INTEGER, INTENT(INOUT) :: N2, N3
      INTEGER :: ND1, ND2, J

      IF (C2.LT.C1) GOTO 991
      IF (C3.LT.C2) GOTO 992

C  IS THERE A 3RD SECTION ?

      IF (CT.GT.C3.AND.N3.GT.1.AND.N3.LT.NT) THEN
C YES. KNOTS IN THIS SECTION:
        ND2=NT-N3+1
      ELSE
C NO.
        ND2=1
        CT=C3
        N3=NT
      ENDIF
C
C  IS THERE A 2ND SECTION
C
      IF (C3.GT.C2.AND.N2.GT.1.AND.N2.LT.N3) THEN
C YES. KNOTS IN THIS SECTION:
        ND1=N3-N2+1
      ELSE
C NO.
        ND1=1
        C3=C2
        N2=N3
      ENDIF

      DO 101 J=1,N2
        C(J)=C1+DBLE(J-1)/DBLE(N2-1)*(C2-C1)
101   CONTINUE
      DO 102 J=N2+1,N3
        C(J)=C2+DBLE(J-N2)/DBLE(N3-N2)*(C3-C2)
102   CONTINUE
      DO 103 J=N3+1,NT
        C(J)=C3+DBLE(J-N3)/DBLE(NT-N3)*(CT-C3)
103   CONTINUE
C
      RETURN
C
991   CONTINUE
      IF (IL.EQ.1) THEN
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: 1ST GRID.  RGA > RIA ?'
        WRITE (iunout,*) 'RI1,RGA = ',C1,C2
      ELSEIF (IL.EQ.2) THEN
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: 2ND GRID.  YGA > YIA ?'
        WRITE (iunout,*) 'YIA,YGA = ',C1,C2
      ELSEIF (IL.EQ.3) THEN
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: 3RD GRID.  ZGA > ZIA ?'
        WRITE (iunout,*) 'ZIA,ZGA = ',C1,C2
      ELSE
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: IN GRID_1   C2 > C1  ?'
        WRITE (iunout,*) 'C1,C2 = ',C1,C2
      ENDIF
      CALL EXIT_OWN(1)
992   CONTINUE
      IF (IL.EQ.1) THEN
        WRITE (iunout,*) 
     .   'GRID DATA INCONSISTENCY: 1ST GRID.  RAA > RGA ?'
        WRITE (iunout,*) 'RGA,RAA = ',C2,C3
      ELSEIF (IL.EQ.2) THEN
        WRITE (iunout,*) 
     .   'GRID DATA INCONSISTENCY: 2ND GRID.  YAA > YGA ?'
        WRITE (iunout,*) 'YGA,YAA = ',C2,C3
      ELSEIF (IL.EQ.3) THEN
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: 3RD GRID.  ZAA > ZGA ?'
        WRITE (iunout,*) 'ZGA,ZAA = ',C2,C3
      ELSE
        WRITE (iunout,*) 
     .    'GRID DATA INCONSISTENCY: IN GRID_1   C3 > C2  ?'
        WRITE (iunout,*) 'C2,C3 = ',C2,C3
      ENDIF
      CALL EXIT_OWN(1)
      RETURN
      END
C ===== SOURCE: grid.f
C
!pb  3.12.06: allow INDGRD /= 6 for NLTET option
!pb  3.12.06: specify NGITT in case of NLTET
!pb           initialize XDIFF=0
      SUBROUTINE GRID (IND)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CTRCEI
      USE CGEOM
      USE CTETRA
      USE CLGIN
      USE CTRIG

      IMPLICIT NONE

      REAL(DP) :: PC1(3), EDGELEN(6)
      REAL(DP) :: ELPARM, X1, X2, SY, Y1, Y2, SX, AELL, X3, Y3, X4, Y4,
     .          VPXX, RN, RRN, FN, VPYY, PLABS, XNORM, VPX, VPY, QUOTI,
     .          GESFL, FRING, CONST, RRR, FL, FR, RL, RR, RRL, XD,
     .          PLEN, XDIFF, RORIG, XS3, PLABS2, PLABS3, XD1, YD,
     .          XS, PLABS1, YD1, XS2, XD3, YD3, XS1, XD2, YD2, R, PIN,
     .          POUT, EX1, SDSD, XX1, XX2, YY1, YY2, DSD, COM, S, SQ
      REAL(DP), EXTERNAL :: ARTRI3
      INTEGER :: ITSIDE(3,4)
      INTEGER :: IP, IRP, IPP, IT, KDN, KUP, NCELL, IR, I, K, IUP, IDN,
     .           J, IND, NLOCAL, ND, IC3, IC4, ITET, IC1, IC2, IC, NLJ,
     .           IFLAG, ISTS, IECKE2, MSURFG, IS, IT1, NCELL1, NSRFTR
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/
C STATEMENT FUNCTION FOR GRID PARAMETERS FOR LEVGEO=2 OPTION
      ELPARM(R,PIN,POUT,EX1)=(PIN-POUT)*(1.-R**EX1)**1.+POUT
C
      GOTO(100,200,300),IND
C
C   RADIAL GRID
C
100   CONTINUE
C
      IF (NR1ST.LT.2) RETURN
C
      IF (LEVGEO.EQ.1) THEN
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.1
C
        IF (INDGRD(IND).LE.4) THEN
C** USE ONE OF THE EIRENE DEFAULT GRID OPTIONS
          CALL GRID_1(RSURF,NR1ST,NRSEP,NRPLG,RIA,RGA,RAA,RRA,1)
        ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE RADIAL GRID DATA FROM USER SUPPLIED SUBROUTINE
C         CALL PROUSR (RSURF,2+4*NPLS+3,0._DP,0._DP,0._DP,0._DP,
C    .                 0._DP,0._DP,0._DP,NR1ST)
        ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE RADIAL GRID DATA FROM INTERFACE
C         CALL PROFR (RSURF,6+5*NPLS+NAIN,1,1,NR1ST)
        ENDIF
C  SET DERIVED GRID DATA FOR LEVGEO = 1 OPTION
C
C  NOTHING TO BE DONE HERE
C
        IF (TRCGRD) THEN
          CALL LEER(1)
          WRITE (iunout,*) 'GRIDPOINTS IN X DIRECTION '
          CALL LEER(1)
          CALL MASRR1('  N, RSURF ',RSURF,NR1ST,3)
          CALL LEER(2)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.2
C
        IF (INDGRD(IND).LE.4) THEN
C** USE ONE OF THE EIRENE DEFAULT GRID OPTIONS
          IF (RRA.GT.RAA) THEN
            NLOCAL=NR1STM
            RSURF(NR1ST)=RRA
            ELL(NR1ST)=ELLCH
            EP1(NR1ST)=EP1CH
            TRI(NR1ST)=TRICH
          ELSE
            NLOCAL=NR1ST
          ENDIF
C
          IF (INDGRD(IND).EQ.1) THEN
            ND=NRSEP
            DO 105 J=1,ND
              RSURF(J)=RIA+DBLE(J-1)/DBLE(ND-1)*(RGA-RIA)
105         CONTINUE
            DO 106 J=ND+1,NLOCAL
              RSURF(J)=RGA+DBLE(J-ND)/DBLE(NLOCAL-ND)*(RAA-RGA)
106         CONTINUE
          ELSEIF (INDGRD(IND).EQ.2) THEN
C** RADIAL GRID WITH CONSTANT AREA
            IF (NLCRC) THEN
              GESFL=(RAA*RAA-RIA*RIA)*PIA
              FRING=GESFL/(NLOCAL-1)
              RSURF(1)=RIA
              RSURF(NLOCAL)=RAA
              DO 103 J=2,NLOCAL-1
                RSURF(J)=RIA+SQRT((J-1)*FRING*PIAI)
103           CONTINUE
            ELSEIF (NLELL) THEN
              GESFL=(RAA*RAA*ELLOT-RIA*RIA*ELLIN)*PIA
              FRING=GESFL/(NLOCAL-1)
              RSURF(1)=RIA
              RSURF(NLOCAL)=RAA
              DO 104 J=2,NLOCAL-1
C  SOLVE RSURF**2*ELL(RSURF)-(J-1)*FRING/PIA=0., RSURF(J-1)<RSURF<RAA
                CONST=(J-1)*FRING*PIAI
                RL=RSURF(J-1)+EPS30
                RR=RAA
108             RRL=(RL-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
109             RRR=(RR-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
                FL=RL**2*ELPARM(RRL,ELLIN,ELLOT,EXELL)-CONST
                FR=RR**2*ELPARM(RRR,ELLIN,ELLOT,EXELL)-CONST
                QUOTI=(RR-RL)/(FR-FL)
                RN=-FL*QUOTI+RL
                RRN=(RN-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
                FN=RN**2*ELPARM(RRN,ELLIN,ELLOT,EXELL)-CONST
                IF (ABS(FN/CONST).LT.EPS6) THEN
                  RSURF(J)=RN
                  GOTO 104
                ELSEIF (FN.LT.0.D0) THEN
                  RL=RN
                  GOTO 108
                ELSEIF (FN.GT.0.D0) THEN
                  RR=RN
                  GOTO 109
                ENDIF
104           CONTINUE
            ELSEIF (NLTRI) THEN
              GESFL=((RAA*RAA-2.*TRIOT*TRIOT)*ELLOT-
     .               (RIA*RIA-2.*TRIIN*TRIIN)*ELLIN)*PIA
              FRING=GESFL/(NLOCAL-1)
              RSURF(1)=RIA
              RSURF(NLOCAL)=RAA
              WRITE (iunout,*) 'NLTRI NOT READY IN GRID'
              CALL EXIT_OWN(1)
            ENDIF
          ELSE
            WRITE (iunout,*) 'INVALID OPTION ENCOUNTERED IN SUBR. GRID'
            WRITE (iunout,*) 'INDGRD(IND),NLCRC ',INDGRD(IND),NLCRC
            WRITE (iunout,*) 'EXIT CALLED FROM SUBR. GRID'
            CALL EXIT_OWN(1)
          ENDIF
        ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE RADIAL GRID DATA FROM USER SUPPLIED SUBROUTINE
          CALL PROUSR (RSURF,2+4*NPLS+3,0._DP,0._DP,0._DP,0._DP,
     .                 0._DP,0._DP,0._DP,NR1ST)
        ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE RADIAL GRID DATA FROM INTERFACE
C         CALL PROFR (RSURF,2+4*NPLS+3,1,1,NR1ST)
        ENDIF
C
C  SET DATA FOR GRID OF ELLIPSES
C
        IF (INDGRD(IND).LE.2) THEN
          IF (NLCRC) THEN
            DO 111 J=1,NR1ST
              EP1(J)=0.
              ELL(J)=1.
              TRI(J)=0.
111         CONTINUE
          ELSEIF (NLELL) THEN
            EP1(1)=EP1IN
            ELL(1)=ELLIN
            DO 112 J=2,NLOCAL
              RR=(RSURF(J)-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
              EP1(J)=ELPARM(RR,EP1IN,EP1OT,EXEP1)
              ELL(J)=ELPARM(RR,ELLIN,ELLOT,EXELL)
112         CONTINUE
            DO 113 J=1,NR1ST
              TRI(J)=0.
113         CONTINUE
          ELSEIF (NLTRI) THEN
            EP1(1)=EP1IN
            ELL(1)=ELLIN
            TRI(1)=TRIIN
            IF (ABS(TRI(1)/(RSURF(1)+EPS30)).GT.EPS30) THEN
              WRITE (iunout,*) 'FROM SUBR. GRID: '
              WRITE (iunout,*) 'ERROR IN TRIANGULARITY PARAMETERS '
              CALL EXIT_OWN(1)
            ENDIF
            DO 114 J=2,NLOCAL
              RR=(RSURF(J)-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
              EP1(J)=ELPARM(RR,EP1IN,EP1OT,EXEP1)
              ELL(J)=ELPARM(RR,ELLIN,ELLOT,EXELL)
              TRI(J)=ELPARM(RR,TRIIN,TRIOT,EXTRI)
114         CONTINUE
          ENDIF
        ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE ELLIP. GRID DATA FROM USER SUPPLIED SUBROUTINE
C         CALL PROUSR (EP1,2+4*NPLS+?,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,
C         CALL PROUSR (ELL,2+4*NPLS+?,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,
C         CALL PROUSR (TRI,2+4*NPLS+?,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,
        ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE ELLIP. GRID DATA FROM INTERFACE
C         CALL PROFR (EP1,2+4*NPLS+?,1,1,NR1ST)
C         CALL PROFR (ELL,2+4*NPLS+?,1,1,NR1ST)
C         CALL PROFR (TRI,2+4*NPLS+?,1,1,NR1ST)
        ENDIF
C
C  SET DERIVED GRID DATA FOR LEVGEO = 2 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C
        DO 115 J=1,NR1ST
          RQ(J)=RSURF(J)*RSURF(J)
          ELLQ(J)=ELL(J)*ELL(J)
115     CONTINUE
C
        IF (TRCGRD) THEN
          CALL MASRR4('  N, RSURF,EP1,ELL,TRI',
     .                      RSURF,EP1,ELL,TRI,NR1ST)
          CALL LEER(2)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.3) THEN
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.3
C
        IF (INDGRD(IND).LE.4) THEN
C  ALL POLYGON DATA HAVE BEEN READ FROM INPUT FILE, NOTHING ELSE
C  TO BE DONE HERE
        ELSEIF (INDGRD(IND).EQ.5) THEN
C*** POLYGON DATA NOT JET AVAILABLE FROM PROUSR (INDGRD.EQ.5 OPTION)
        ELSEIF (INDGRD(IND).EQ.6) THEN
C*** GEOMETRICAL DATA ARE SET IN IF0COP NOTHING TO BE DONE HERE
        ENDIF
C
        IF (PLREFL.GT.0._DP) THEN
C  SET ONE ADDITIONAL OUTERMOST POLYGON, DISTANCE PLREFL (CM) FROM THE
C  OUTERMOST POLYGON SPECIFIED SO FAR
          I=NR1STM
          DO 133 J=1,NPPLG
            DO 132 K=NPOINT(1,J),NPOINT(2,J)
              VPXX=XPOL(I,K)-XPOL(I-1,K)
              VPYY=YPOL(I,K)-YPOL(I-1,K)
              XNORM=SQRT(VPXX**2+VPYY**2)
              VPX=VPXX/XNORM*PLREFL
              VPY=VPYY/XNORM*PLREFL
              XPOL(I+1,K)=XPOL(I,K)+VPX
              YPOL(I+1,K)=YPOL(I,K)+VPY
132         CONTINUE
133       CONTINUE
        ENDIF
C
        IF (XPCOR.NE.0._DP) THEN
C  SHIFT WHOLE POLYGON MESH IN X DIRECTION BY XPCOR (CM)
C
          DO 135 I=1,NR1ST
            DO 135 J=1,NPPLG
              DO 135 K=NPOINT(1,J),NPOINT(2,J)
                XPOL(I,K)=XPOL(I,K)+XPCOR
135       CONTINUE
        ENDIF
C
        IF (YPCOR.NE.0._DP) THEN
C  SHIFT WHOLE POLYGON MESH IN Y DIRECTION BY YPCOR (CM)
C
          DO 136 I=1,NR1ST
            DO 136 J=1,NPPLG
              DO 136 K=NPOINT(1,J),NPOINT(2,J)
                YPOL(I,K)=YPOL(I,K)+YPCOR
136       CONTINUE
        ENDIF
C
C  SET DERIVED GRID DATA FOR LEVGEO = 3 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C
        DO 140 I=1,NR1ST
          DO 140 J=1,NPPLG
            DO 140 K=NPOINT(1,J),NPOINT(2,J)-1
              VPLX(I,K)=XPOL(I,K+1)-XPOL(I,K)
              VPLY(I,K)=YPOL(I,K+1)-YPOL(I,K)
140     CONTINUE
        DO 141 I=1,NR1ST
          DO 141 J=1,NPPLG
            IF (J.EQ.1) THEN
              DO 142 K=1,NPOINT(1,1)
                BGL(I,K)=0.
142           CONTINUE
            ELSE
              DO 143 K=NPOINT(2,J-1),NPOINT(1,J)
                BGL(I,K)=BGL(I,NPOINT(2,J-1))
143           CONTINUE
            ENDIF
            DO 141 K=NPOINT(1,J)+1,NPOINT(2,J)
              BGL(I,K)=BGL(I,K-1)+SQRT(VPLX(I,K-1)**2+VPLY(I,K-1)**2)
141     CONTINUE
C
C   CALCULATE THE OUTER NORMALS OF POLYGONS
C
        DO 144 I=1,NR1ST
          DO 144 J=1,NPPLG
            DO 144 K=NPOINT(1,J),NPOINT(2,J)-1
              PLABS=SQRT(VPLX(I,K)**2+VPLY(I,K)**2)
              PLNX(I,K)=VPLY(I,K)/(PLABS+EPS60)
              PLNY(I,K)=-VPLX(I,K)/(PLABS+EPS60)
144     CONTINUE
C
        DO 147 I=1,NR1ST
          IUP=I+1
          IDN=I
          IF (IUP.GT.NR1ST) THEN
            IUP=I
            IDN=I-1
          ENDIF
          DO 147 J=1,NPPLG
          DO 147 K=NPOINT(1,J),NPOINT(2,J)-1
146         XD=XPOL(IUP,K+1)-XPOL(IDN,K+1)
            YD=YPOL(IUP,K+1)-YPOL(IDN,K+1)
            IF (XD*XD+YD*YD.LT.EPS30) THEN
              IF (IUP.LT.NR1ST) THEN
                IUP=IUP+1
              ELSE
                IDN=IDN-1
              ENDIF
              GOTO 146
            ENDIF
            XS=SIGN(1._DP,XD*PLNX(I,K)+YD*PLNY(I,K))
            PLNX(I,K)=PLNX(I,K)*XS
            PLNY(I,K)=PLNY(I,K)*XS
147     CONTINUE
C
        IF (TRCGRD) THEN
          WRITE (iunout,*) ' NO. OF VALID PARTS = ',NPPLG
          DO 155 J=1,NR1ST
            WRITE (iunout,*) ' POLYGON NO. J = ',J
            DO 156 K=1,NPPLG
              WRITE (iunout,*) 'IA = ',NPOINT(1,K),' IE = ',NPOINT(2,K)
              WRITE (iunout,'(/1X,1P,6E12.4)') (XPOL(J,I),YPOL(J,I),
     .                                   I=NPOINT(1,K),NPOINT(2,K))
156         CONTINUE
155       CONTINUE
          CALL LEER(2)
          WRITE (iunout,*) 
     .      'ARCLENGTH BGL(I,K) OF RADIAL SURFACES AT Z=0.'
          DO 153 I=1,NR1ST
            WRITE (iunout,*) 'I = ',I
            WRITE (iunout,'(/1X,1P,6E12.4)') (BGL(I,K),K=1,NP2ND)
            CALL LEER(1)
153       CONTINUE
        ENDIF
C
        CALL SNEIGH
C
      ELSEIF (LEVGEO.EQ.4) THEN
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.4
C
!pb        IF (INDGRD(1).NE.6) THEN
!pb          WRITE (iunout,*) 
!PB     .      ' WRONG GRID OPTION SPECIFIED FOR FEM GRID '
!pb          CALL EXIT_OWN(1)
!pb        ENDIF
C
C  SET DERIVED GRID DATA FOR LEVGEO = 4 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C
        DO 165 I=1,NTRII
          VTRIX(1,I)=XTRIAN(NECKE(2,I))-XTRIAN(NECKE(1,I))
          VTRIY(1,I)=YTRIAN(NECKE(2,I))-YTRIAN(NECKE(1,I))
          VTRIX(2,I)=XTRIAN(NECKE(3,I))-XTRIAN(NECKE(2,I))
          VTRIY(2,I)=YTRIAN(NECKE(3,I))-YTRIAN(NECKE(2,I))
          VTRIX(3,I)=XTRIAN(NECKE(1,I))-XTRIAN(NECKE(3,I))
          VTRIY(3,I)=YTRIAN(NECKE(1,I))-YTRIAN(NECKE(3,I))
165     CONTINUE
C
C
C   CALCULATE THE OUTER NORMALS OF TRIANGLES
C
        DO 161 I=1,NTRII
          PLABS1=SQRT(VTRIX(1,I)**2+VTRIY(1,I)**2)
          PLABS2=SQRT(VTRIX(2,I)**2+VTRIY(2,I)**2)
          PLABS3=SQRT(VTRIX(3,I)**2+VTRIY(3,I)**2)
          PTRIX(1,I)=VTRIY(1,I)/(PLABS1+EPS60)
          PTRIX(2,I)=VTRIY(2,I)/(PLABS2+EPS60)
          PTRIX(3,I)=VTRIY(3,I)/(PLABS3+EPS60)
          PTRIY(1,I)=-VTRIX(1,I)/(PLABS1+EPS60)
          PTRIY(2,I)=-VTRIX(2,I)/(PLABS2+EPS60)
          PTRIY(3,I)=-VTRIX(3,I)/(PLABS3+EPS60)
161     CONTINUE
C
C PTRIX/PTRIY POINT OUT OF TRIANGE FOR A MATHEMATICAL POSITVE
C             ORIENTATION OF TRIANGLE (1-2-3-1: COUNTER CLOCKWISE,
C                                               AS IT MUST BE)
C PTRIX/PTRIY POINT INTO TRIANGE FOR A MATHEMATICAL NEGATIVE
C             ORIENTATION OF TRIANGLE (1-2-3-1: CLOCKWISE,
C                                               ERROR MESSAGE FROM LEARC1)
C
        DO 162 I=1,NTRII
          XD1=VTRIX(1,I)+PTRIX(1,I)
          YD1=VTRIY(1,I)+PTRIY(1,I)
          XS1=SIGN(1._DP,XD1*VTRIX(1,I)+YD1*VTRIY(1,I))
          PTRIX(1,I)=PTRIX(1,I)*XS1
          PTRIY(1,I)=PTRIY(1,I)*XS1
          XD2=VTRIX(2,I)+PTRIX(2,I)
          YD2=VTRIY(2,I)+PTRIY(2,I)
          XS2=SIGN(1._DP,XD2*VTRIX(2,I)+YD2*VTRIY(2,I))
          PTRIX(2,I)=PTRIX(2,I)*XS2
          PTRIY(2,I)=PTRIY(2,I)*XS2
          XD3=VTRIX(3,I)+PTRIX(3,I)
          YD3=VTRIY(3,I)+PTRIY(3,I)
          XS3=SIGN(1._DP,XD3*VTRIX(3,I)+YD3*VTRIY(3,I))
          PTRIX(3,I)=PTRIX(3,I)*XS3
          PTRIY(3,I)=PTRIY(3,I)*XS3
162     CONTINUE
C
C  DETERMINE NGITT FROM NUMBER OF NONDEFAULT STANDARD SURFACES
C
        NGITT = COUNT(INMTI(1:3,1:NTRII) .NE. 0) + 1

        IF (MAXVAL(INMTI(1:3,1:NTRII)) > NSTSI+NLIM) THEN
          WRITE (iunout,*) 
     .      ' WRONG NUMBER OF REFLECTION MODEL SPECIFIED '
          WRITE (iunout,*) ' CHECK DEFINITION OF TRIANGLES ',
     .                'AND THEIR NEIGHBORS '
          CALL EXIT_OWN(1)
        END IF

        DO I = 1, NLIMPS
          NSRFTR = COUNT(INMTI(1:3,1:NTRII) .EQ. I)
          IF (NSRFTR > 0) THEN
            ALLOCATE (SURF_TRIAN(I)%ITRIAS(NSRFTR))
            ALLOCATE (SURF_TRIAN(I)%ITRISI(NSRFTR))
            ALLOCATE (SURF_TRIAN(I)%BGLT(NSRFTR+1))
            SURF_TRIAN(I)%BGLT(1) = 0._DP
          END IF
        END DO

        DO I=1,NTRII
          DO IS = 1, 3
            J = INMTI(IS,I)
            IF ( J .NE. 0) THEN
              SURF_TRIAN(J)%NUMTR = SURF_TRIAN(J)%NUMTR + 1
              SURF_TRIAN(J)%ITRIAS(SURF_TRIAN(J)%NUMTR) = I
              SURF_TRIAN(J)%ITRISI(SURF_TRIAN(J)%NUMTR) = IS
              SURF_TRIAN(J)%BGLT(SURF_TRIAN(J)%NUMTR+1) = 
     .          SURF_TRIAN(J)%BGLT(SURF_TRIAN(J)%NUMTR) + 
     .          SQRT(VTRIX(IS,I)**2+VTRIY(IS,I)**2)
            END IF
          END DO
        END DO

C
        IF (TRCGRD) THEN
          WRITE (iunout,*) ' NUMBER OF TRIANGLES = ',NTRII
          WRITE (iunout,*) ' I,(XTRIAN(J),YTRIAN(J),J=1,3) '
          DO 163 I=1,NTRII
            WRITE (iunout,'(/1X,I4,1X,1P,6E12.4)')
     .                               I,(XTRIAN(NECKE(J,I)),
     .                                  YTRIAN(NECKE(J,I)),J=1,3)
163       CONTINUE
          CALL LEER(2)
          WRITE (iunout,*) ' NGITT SET TO ',NGITT

          DO I=1, NLIMPS
            WRITE (IUNOUT,*) 
            WRITE (IUNOUT,*) ' SURFACE NO. ',I
            WRITE (IUNOUT,'(3A6,A12)') 'J','ITRI','ISIDE','BLGT'
            DO J=1, SURF_TRIAN(I)%NUMTR
              WRITE (IUNOUT,'(3I6,ES12.4)') J, SURF_TRIAN(I)%ITRIAS(J), 
     .                                         SURF_TRIAN(I)%ITRISI(J),
     .                                         SURF_TRIAN(I)%BGLT(J+1)
            END DO
          END DO
        ENDIF
C
C
      ELSEIF (LEVGEO.EQ.5) THEN
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.5
C
!pb        IF (INDGRD(1).NE.6) THEN
!pb          WRITE (iunout,*) ' WRONG GRID OPTION SPECIFIED FOR',
!pb     .                ' TETRAHEDRON GRID '
!pb          CALL EXIT_OWN(1)
!pb        ENDIF
C  GRID DATA FOR TETRAHEDRONS ARE SET IN COUPLING ROUTINE
C  NOTHING TO BE DONE HERE
C
C  SET DERIVED GRID DATA FOR LEVGEO = 6 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C

        DO ITET=1,NTET
          IC1 = NTECK(1,ITET)
          IC2 = NTECK(2,ITET)
          IC3 = NTECK(3,ITET)
          IC4 = NTECK(4,ITET)
          RINCRC(1:4,ITET) = 1._DP
C  CALCULATE DIRECTIONS OF EDGES
C  EDGE  1-2
          VTETX(1,ITET) = XTETRA(IC2) - XTETRA(IC1)
          VTETY(1,ITET) = YTETRA(IC2) - YTETRA(IC1)
          VTETZ(1,ITET) = ZTETRA(IC2) - ZTETRA(IC1)
C  EDGE  2-3
          VTETX(2,ITET) = XTETRA(IC3) - XTETRA(IC2)
          VTETY(2,ITET) = YTETRA(IC3) - YTETRA(IC2)
          VTETZ(2,ITET) = ZTETRA(IC3) - ZTETRA(IC2)
C  EDGE  3-1
          VTETX(3,ITET) = XTETRA(IC1) - XTETRA(IC3)
          VTETY(3,ITET) = YTETRA(IC1) - YTETRA(IC3)
          VTETZ(3,ITET) = ZTETRA(IC1) - ZTETRA(IC3)
C  EDGE  1-4
          VTETX(4,ITET) = XTETRA(IC4) - XTETRA(IC1)
          VTETY(4,ITET) = YTETRA(IC4) - YTETRA(IC1)
          VTETZ(4,ITET) = ZTETRA(IC4) - ZTETRA(IC1)
C  EDGE  2-4
          VTETX(5,ITET) = XTETRA(IC4) - XTETRA(IC2)
          VTETY(5,ITET) = YTETRA(IC4) - YTETRA(IC2)
          VTETZ(5,ITET) = ZTETRA(IC4) - ZTETRA(IC2)
C  EDGE  3-4
          VTETX(6,ITET) = XTETRA(IC4) - XTETRA(IC3)
          VTETY(6,ITET) = YTETRA(IC4) - YTETRA(IC3)
          VTETZ(6,ITET) = ZTETRA(IC4) - ZTETRA(IC3)

          EDGELEN(1:6) = SQRT(VTETX(1:6,ITET)**2 +
     .                        VTETY(1:6,ITET)**2 +
     .                        VTETZ(1:6,ITET)**2)
C  CALCULATE THE OUTER NORMALS OF TETRAHEDRONS
C  SIDE 1-2-3
          PTETX(1,ITET) = VTETY(3,ITET)*VTETZ(1,ITET) -
     .                    VTETZ(3,ITET)*VTETY(1,ITET)
          PTETY(1,ITET) = VTETZ(3,ITET)*VTETX(1,ITET) -
     .                    VTETX(3,ITET)*VTETZ(1,ITET)
          PTETZ(1,ITET) = VTETX(3,ITET)*VTETY(1,ITET) -
     .                    VTETY(3,ITET)*VTETX(1,ITET)
          S = 0.5_DP *( EDGELEN(1) + EDGELEN(2) + EDGELEN(3) )
          SQ = SQRT( (S-EDGELEN(1)) *
     .               (S-EDGELEN(2)) * (S-EDGELEN(3)) / S)
          IF (SQ > EPS30) RINCRC(1,ITET) = 1._DP / SQ
C  SIDE 1-2-4
          PTETX(2,ITET) = VTETY(4,ITET)*VTETZ(1,ITET) -
     .                    VTETZ(4,ITET)*VTETY(1,ITET)
          PTETY(2,ITET) = VTETZ(4,ITET)*VTETX(1,ITET) -
     .                    VTETX(4,ITET)*VTETZ(1,ITET)
          PTETZ(2,ITET) = VTETX(4,ITET)*VTETY(1,ITET) -
     .                    VTETY(4,ITET)*VTETX(1,ITET)
          S =  0.5_DP *( EDGELEN(1) + EDGELEN(5) + EDGELEN(4) )
          SQ = SQRT( (S-EDGELEN(1)) *
     .               (S-EDGELEN(5)) * (S-EDGELEN(4)) / S)
          IF (SQ > EPS30) RINCRC(2,ITET) = 1._DP / SQ
C  SIDE 2-3-4
          PTETX(3,ITET) = VTETY(5,ITET)*VTETZ(2,ITET) -
     .                    VTETZ(5,ITET)*VTETY(2,ITET)
          PTETY(3,ITET) = VTETZ(5,ITET)*VTETX(2,ITET) -
     .                    VTETX(5,ITET)*VTETZ(2,ITET)
          PTETZ(3,ITET) = VTETX(5,ITET)*VTETY(2,ITET) -
     .                    VTETY(5,ITET)*VTETX(2,ITET)
          S =  0.5_DP *( EDGELEN(2) + EDGELEN(6) + EDGELEN(5) )
          SQ = SQRT( (S-EDGELEN(2)) *
     .               (S-EDGELEN(6)) * (S-EDGELEN(5)) / S)
          IF (SQ > EPS30) RINCRC(3,ITET) = 1._DP / SQ
C  SIDE 3-1-4
          PTETX(4,ITET) = VTETY(4,ITET)*VTETZ(3,ITET) -
     .                    VTETZ(4,ITET)*VTETY(3,ITET)
          PTETY(4,ITET) = VTETZ(4,ITET)*VTETX(3,ITET) -
     .                    VTETX(4,ITET)*VTETZ(3,ITET)
          PTETZ(4,ITET) = VTETX(4,ITET)*VTETY(3,ITET) -
     .                    VTETY(4,ITET)*VTETX(3,ITET)
          S =  0.5_DP *( EDGELEN(3) + EDGELEN(4) + EDGELEN(6) )
          SQ = SQRT( (S-EDGELEN(3)) *
     .               (S-EDGELEN(4)) * (S-EDGELEN(6)) / S)
          IF (SQ > EPS30) RINCRC(4,ITET) = 1._DP / SQ

          DO J=1,4
            PLEN=SQRT(PTETX(J,ITET)**2+PTETY(J,ITET)**2+
     .                PTETZ(J,ITET)**2)+EPS60
            PTETX(J,ITET)=PTETX(J,ITET)/PLEN
            PTETY(J,ITET)=PTETY(J,ITET)/PLEN
            PTETZ(J,ITET)=PTETZ(J,ITET)/PLEN
          END DO
        END DO

        NGITT = COUNT(INMTIT(1:4,1:NTET) .NE. 0) + 1

        CALL SUCHE_NACHBARN

        IC=0
        NTET_COLLAPS=0
        DO ITET=1,NTET
          DO IS=1,4
            IF ((NTBAR(IS,ITET) == 0) .AND. (INMTIT(IS,ITET) == 0)) THEN
              IC=IC+1
              WRITE (iunout,*) ' TETRAHEDRON WITH NO NEIGHBORS AND NO ',
     .                    'REFLECTION MODEL FOUND '
              WRITE (iunout,*) ' ITET = ',ITET,' ISIDE = ',IS
              write (iunout,*) nteck(itside(1,is),itet),
     .                    nteck(itside(2,is),itet),
     .                    nteck(itside(3,is),itet)
            END IF
            IF (NTBAR(IS,ITET) < 0) THEN
              IF (SUM(NTBAR(1:4,ITET)) > -4) THEN
                WRITE (iunout,*) 
     .            ' TETRAHEDRON WITH NEIGHBOR -1 DETECTED '
                WRITE (iunout,*) ' ITET = ',ITET,' ISIDE = ',IS
                WRITE (iunout,*) ' NTBAR(ITET) = ',NTBAR(1:4,ITET)
                IC=IC+1
              ELSE
                NTET_COLLAPS=NTET_COLLAPS+1
!pb                WRITE(iunout,*) ' COLLAPSED TETRAHEDRON ITET = ',ITET
                EXIT
              END IF
            END IF
          END DO
        END DO

        IF (IC > 0) CALL EXIT_OWN(1)
C
        IF (TRCGRD) THEN
          WRITE (iunout,*) ' NUMBER OF COORDINATES = ',NCOOR
          WRITE (iunout,*) ' I,(XETRA(J),YTETRA(J),ZTETRA(J),J=1,3) '
          DO I=1,NCOOR
            WRITE (iunout,'(1X,I4,1X,6ES12.4)')
     .             I,XTETRA(I),YTETRA(I),ZTETRA(I)
          END DO
          CALL LEER(2)

          WRITE (iunout,*) ' NUMBER OF TETRAHEDRONS = ',NTET
          DO ITET=1,NTET
            WRITE (iunout,*)
            WRITE (iunout,*) ' TETRAEDER ',ITET
            DO J=1,4
              IC=NTECK(J,ITET)
              WRITE (iunout,'(1X,I6,3ES12.4)')
     .               IC,XTETRA(IC),YTETRA(IC),ZTETRA(IC)
            END DO
            DO J=1,4
              WRITE (iunout,'(A,I3,A,3I6,4X,A,I6,A,I6)')
     .              ' SIDE ',J,': ',
     .                NTECK(ITSIDE(1,J),ITET),
     .                NTECK(ITSIDE(2,J),ITET),
     .                NTECK(ITSIDE(3,J),ITET),
     .              ' NEIGHBOR ',NTBAR(J,ITET),
     .              ' SIDE ',NTSEITE(J,ITET)
            END DO
            WRITE (iunout,*) 'OUTER NORMALS '
            DO J=1,4
              WRITE (iunout,'(7X,3ES12.4)')
     .              PTETX(J,ITET),PTETY(J,ITET),PTETZ(J,ITET)
            END DO
          END DO

C  CHECK OUTER NORMALS
          DO ITET = 1,NTET
            DO J=1,4
             IF (NTBAR(J,ITET) /= 0) THEN
               PC1(1:3) = (/
     .         ABS(PTETX(J,ITET)+PTETX(NTSEITE(J,ITET),NTBAR(J,ITET))),
     .         ABS(PTETY(J,ITET)+PTETY(NTSEITE(J,ITET),NTBAR(J,ITET))),
     .         ABS(PTETZ(J,ITET)+PTETZ(NTSEITE(J,ITET),NTBAR(J,ITET)))/)
               IF (ANY(PC1 > 1.E-6)) THEN
                WRITE(iunout,*) ' PROBLEM WITH TETRAHEDRON ',ITET,
     .                          ' SIDE ',J
                WRITE(iunout,*) ' NORMALS DO NOT MATCH '
                WRITE(iunout,*) PTETX(J,ITET),PTETY(J,ITET),
     .                          PTETZ(J,ITET)
                WRITE(iunout,*) PTETX(NTSEITE(J,ITET),NTBAR(J,ITET)),
     .                     PTETY(NTSEITE(J,ITET),NTBAR(J,ITET)),
     .                     PTETZ(NTSEITE(J,ITET),NTBAR(J,ITET))
               END IF
             END IF
            END DO
          END DO

        ENDIF
C
      ELSEIF (LEVGEO.EQ.6) THEN
C
C  GENERAL GEOMETRY OPTION: NOTHING TO DONE HERE
C
      ENDIF
C
C  SET GEOMETRICAL CONTANTS FOR IGNORABLE Y OR POLOIDAL CO-ORDINATE
C  THESE MAY BE REVISED IF A 2ND (Y- OR POL.) GRID IS DEFINED BELOW
C
      IF (LEVGEO.EQ.1) THEN
        YDF=YAA-YIA
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.2) THEN
        YDF=(YAA-YIA)*DEGRAD
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.3) THEN
        YDF=1.
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.4) THEN
        YDF=1.
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.5) THEN
        YDF=1.
      ELSEIF (LEVGEO.EQ.6) THEN
        YDF=1.
C
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
C
      ENDIF
      IF (YDF.LE.0._DP) GOTO 992
C
      IF (TRCGRD.AND..NOT.NLPOL) THEN
        CALL LEER(2)
        WRITE (iunout,*) 'CONSTANTS FOR POLOIDAL OR Y DIRECTION'
        CALL MASR3('YDF,YIA,YAA=            ',YDF,YIA,YAA)
        CALL LEER(1)
      ENDIF
C
C  SET GEOMETRICAL CONTANTS FOR IGNORABLE Z CO-ORDINATE
C  THESE MAY BE REVISED IF A 3RD (Z- OR TOR.) GRID IS DEFINED BELOW
C
C  A) IN TOROIDAL APPROXIMATION:
C
      IF (NLTRA) THEN
        IF (NTTRA.LE.3.OR.ROA.LT.0._DP) GOTO 991
        IF (LEVGEO.EQ.1) THEN
          XDIFF=0.
        ELSEIF (LEVGEO.EQ.2) THEN
          XDIFF=EP1OT
        ELSEIF (LEVGEO.EQ.3) THEN
          XDIFF=0.
        ELSEIF (LEVGEO.EQ.4) THEN
          XDIFF=0.
        ELSEIF (LEVGEO.EQ.5) THEN
          XDIFF=0.
C
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
C
        ENDIF
C
C  TOROIDAL ANGLE, IN RADIANS
        ZDF=(ZAA-ZIA)*DEGRAD
C
C  ALPHA: HALF OF THE ANGLE INCREMENT IN EQUIDISTANT TOROIDAL ANGLE GRID
        ALPHA=0.5*(ZDF/DBLE(NTTRAM))
        TANAL=TAN(ALPHA)
        SINAL=SIN(2.*ALPHA)
        COSAL=COS(2.*ALPHA)
C
        DPHI=1./(2.*ALPHA)
C
C  ROA IS THE LARGE RADIUS OF THE TORUS
C  (RMTOR,0,0) IS THE ORIGIN OF LOCAL CO-ORDINATE SYSTEM IN
C              EACH TOROIDAL CELL
C  RMTOR SUCH THAT VOLUME OF TORUS = VOLUME OF THE NTTRAM SEGMENTS
C  AT PRESENT: FULLFILLED AT SURFACE DEFINED BY (RAA,EP1OT,ELLOT)
C              OR AT A POLYGON WITH XDIFF=0. (IF THERE IS ONE)
        RMTOR=(ROA+XDIFF)*ALPHA/TANAL-XDIFF
        ZHALF=ALPHA
        ZFULL=ZHALF*2.
C  SET ZSURF EVEN IF NLTOR=FALSE, FOR 3D GEOMETRY PLOTS
        DO 170 J=1,NTTRA
          ZSURF(J)=ZIA*DEGRAD+(J-1)/DPHI
170     CONTINUE
        DO 172 J=1,NTTRAM
          ZZONE(J)=0.5*(ZSURF(J)+ZSURF(J+1))
172     CONTINUE
        RORIG=RMTOR
C
C  B) IN CYLIND. APPROXIMATION:
C     ROA AND RMTOR ARE IRRELEVANT IN THIS CASE, AND ARE NOT DEFINED
C
      ELSEIF (NLTRZ) THEN
        ZDF=ZAA-ZIA
        IF (ZDF.LE.0._DP) GOTO 991
        ZSURF(1)=ZIA
        ZZONE(1)=(ZAA+ZIA)*0.5
        DPHI=1./ZDF
C
        RORIG=0.
C
C  C) IN TORUS CO-ORDINATES
C
      ELSEIF (NLTRT) THEN
        ZDF=(ZAA-ZIA)*DEGRAD
        IF (ZDF.LE.0._DP) GOTO 991
        ZSURF(1)=ZIA
        ZZONE(1)=(ZAA+ZIA)*0.5
        DPHI=1./ZDF
C
        RORIG=0.
C
      ELSE
        WRITE (iunout,*) ' ERROR IN INPUT DATA! '
        WRITE (iunout,*) ' NLTRA OR NLTRZ OR NLTRT MUST BE .TRUE. '
        CALL EXIT_OWN(1)
      ENDIF
C
      IF (TRCGRD) THEN
        CALL LEER(2)
        IF (.NOT.NLTOR) THEN
          WRITE (iunout,*) 'CONSTANTS FOR TOROIDAL OR Z DIRECTION'
          CALL MASR3('ZDF,ZIA,ZAA=            ',ZDF,ZIA,ZAA)
        ENDIF
        IF (NLTRA) THEN
          WRITE (iunout,*) 'ROA,RMTOR= ',ROA,RMTOR
          IF (.NOT.NLTOR) THEN
            CALL MASRR1 (' N,  ZSURF ',ZSURF,NTTRA,3)
            CALL MASRR1 (' N,  ZZONE ',ZZONE,NTTRAM,3)
          ENDIF
          CALL LEER(2)
        ENDIF
        CALL LEER(1)
      ENDIF
C
C  SET SURFACE AREA OF NON DEFAULT STANDARD SURFACES
C
      IF (LEVGEO.EQ.1) THEN
        DO 180 ISTS=1,NSTSI
          IF (INUMP(ISTS,1).NE.0) THEN
            IR=INUMP(ISTS,1)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=YDF*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=YDF*ZDF*(RSURF(IR)+RMTOR)
            ELSEIF (NLTRT) THEN
            ENDIF
          ENDIF
180     CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
        DO ISTS=1,NSTSI
          NLJ=NLIM+ISTS
          SAREA(NLJ)=0.
        ENDDO
        DO I=1,NTRII
          DO J=1,3
            IECKE2 = J+1
            IF (IECKE2 .EQ. 4) IECKE2 = 1
            ISTS=ABS(INMTI(J,I))
            IF (ISTS .GT. NLIM) THEN
C  SEITE J VON DREIECK I GEHOERT ZUM RAND ISTS
              XX1=XTRIAN(NECKE(J,I))
              YY1=YTRIAN(NECKE(J,I))
              XX2=XTRIAN(NECKE(IECKE2,I))
              YY2=YTRIAN(NECKE(IECKE2,I))
              DSD=((XX1-XX2)**2+(YY1-YY2)**2)**0.5
              IF (NLTRA) THEN
                XX1=XX1+RMTOR
                XX2=XX2+RMTOR
                COM=0.5*(XX1+XX2)
                DSD=DSD*COM*2.*PIA
              ELSE
                DSD=DSD*ZDF
              ENDIF
              IF (NLMPGS.GT.NLIMPS) THEN
                MSURFG=NLIM+NSTS+INSPAT(J,I)
                SAREA(MSURFG)=DSD
              END IF
              SAREA(ISTS)=SAREA(ISTS)+DSD
            ENDIF
          ENDDO
        ENDDO
      ELSEIF (LEVGEO.EQ.5) THEN
        DO ISTS=1,NSTSI
          NLJ=NLIM+ISTS
          SAREA(NLJ)=0.
        ENDDO
        DO ITET=1,NTET
          DO J=1,4
            ISTS=ABS(INMTIT(J,ITET))
            IF (ISTS .GT. 0) THEN
              IC1=NTECK(ITSIDE(1,J),ITET)
              IC2=NTECK(ITSIDE(2,J),ITET)
              IC3=NTECK(ITSIDE(3,J),ITET)
!pb              NLJ=NLIM+ISTS   ! INMTIT CHANGED IN INFCOP
              SAREA(ISTS)=SAREA(ISTS)+
     .                   ARTRI3(XTETRA(IC1),YTETRA(IC1),ZTETRA(IC1),
     .                          XTETRA(IC2),YTETRA(IC2),ZTETRA(IC2),
     .                          XTETRA(IC3),YTETRA(IC3),ZTETRA(IC3))
            END IF
          END DO
        END DO
C     ELSEIF (LEVGEO.EQ....) THEN
      ENDIF

!  SET NSTGRD FOR AVERAGING CELLS

      IR = NR1ST
      DO IT = 1, NT3RD
        DO IP = 1, NP2ND
          NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          NSTGRD(NCELL) = 3
        END DO
      END DO

C
      RETURN
C
C   POLOIDAL OR Y-GRID
C
200   CONTINUE
C
C  IF NLSYMP, Y-GRID MUST BE SYMMETRIC: PSURF(I)=PSURF(NP2ND-I+1)
C
      IF (LEVGEO.EQ.1) THEN
C   Y-GRID
        IF (INDGRD(IND).LE.4) THEN
          CALL GRID_1(PSURF,NP2ND,NPSEP,NPPLA,YIA,YGA,YAA,YYA,2)
C       ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE Y GRID DATA FROM USER SUPPLIED SUBROUTINE
C TO BE WRITTEN
C         CALL PROUSR (PSURF,2+4*NPLS+3,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._D
C       ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE Y GRID DATA FROM INTERFACE
C TO BE WRITTEN
C         CALL PROFR (PSURF,2+4*NPLS+5,1,1,NP2ND)
        ENDIF
C
        DO 210 J=1,NP2NDM
210       PHZONE(J)=(PSURF(J)+PSURF(J+1))/2.
C
        IF (TRCGRD) THEN
          CALL LEER(1)
          WRITE (iunout,*) 'GRIDPOINTS IN Y DIRECTION '
          CALL LEER(1)
          CALL MASRR1('  N, PSURF ',PSURF,NP2ND,3)
          CALL LEER(2)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
        IF (INDGRD(IND).EQ.1) THEN
          ND=NPSEP
          DO 211 J=1,ND
            PSURF(J)=(YIA+DBLE((J-1))/DBLE(ND-1)*(YGA-YIA))*DEGRAD
211       CONTINUE
          DO 212 J=ND+1,NP2ND
            PSURF(J)=(YGA+DBLE(J-ND)/DBLE(NP2ND-ND)*(YAA-YGA))*DEGRAD
212       CONTINUE
        ELSEIF (INDGRD(IND).EQ.2) THEN
          DO 213 J=1,NP2ND
            PSURF(J)=(YIA+(J-1)/DBLE(NP2NDM)*(YAA-YIA))*DEGRAD
213       CONTINUE
        ENDIF
        DO 215 J=1,NP2ND
          COSPH(J)=COS(PSURF(J))
          SINPH(J)=SIN(PSURF(J))
215     CONTINUE
C
        IF (TRCGRD) THEN
          CALL LEER(1)
          WRITE (iunout,*) 'GRIDPOINTS IN POLOIDAL DIRECTION '
          CALL LEER(1)
          CALL MASRR1('  N, PSURF ',PSURF,NP2ND,3)
          CALL LEER(2)
        ENDIF
C
C
        NPPLG=1
        NPOINT(1,1)=1
        NPOINT(2,1)=NP2ND
        IFLAG=2
        DO 1240 IR=1,NR1STM
          IRP=IR+1
          DO 1250 IP=1,NP2NDM
            IPP=IP+1
            CALL ARELLP(EP1(IRP),EP1(IR),ELL(IRP),ELL(IR),
     .                  TRI(IRP),TRI(IR),
     .                  RSURF(IRP),RSURF(IR),PSURF(IPP),PSURF(IP),IFLAG,
     .                  AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
C
            XPOL(IRP,IPP)=X1
            YPOL(IRP,IPP)=Y1
            XPOL(IR,IPP)=X2
            YPOL(IR,IPP)=Y2
            XPOL(IRP,IP)=X3
            YPOL(IRP,IP)=Y3
            XPOL(IR,IP)=X4
            YPOL(IR,IP)=Y4
1250      CONTINUE
1240    CONTINUE
C
        CALL SNEIGH

      ENDIF
C
      IF (LEVGEO.EQ.2.OR.LEVGEO.EQ.3) THEN
C
        IF (TRCGRD) THEN
          DO 219 I=1,NP2ND
            WRITE (iunout,*) ' PERP. POLYGON NO. I = ',I
            WRITE (iunout,*) ' JA = ',1,' JE = ',NR1ST
            WRITE (iunout,'(/1X,1P,6E12.4)') (XPOL(K,I),YPOL(K,I),
     .             K=1,NR1ST)
219       CONTINUE
        ENDIF
C
        DO 220 K=1,NP2ND
          DO 220 I=1,NR1STM
            VVTX(I,K)=XPOL(I+1,K)-XPOL(I,K)
            VVTY(I,K)=YPOL(I+1,K)-YPOL(I,K)
220     CONTINUE
C
        DO 221 K=1,NP2ND
          BGLP(1,K)=0.
          DO 222 I=1,NR1STM
            BGLP(I+1,K)=BGLP(I,K)+SQRT(VVTX(I,K)**2+VVTY(I,K)**2)
222       CONTINUE
221     CONTINUE
C
        IF (TRCGRD) THEN
          CALL LEER(2)
          WRITE (iunout,*) 
     .      'ARCLENGTH BGLP(I,K) OF POLOIDAL SURFACES AT Z=0.'
          DO 223 K=1,NP2ND
            WRITE (iunout,*) 'K = ',K
            WRITE (iunout,'(/1X,1P,6E12.4)') (BGLP(I,K),I=1,NR1ST)
            CALL LEER(1)
223       CONTINUE
        ENDIF
C
C   CALCULATE THE OUTER NORMALS OF POLYGONS
C
        DO 224 K=1,NP2ND
          DO 225 I=1,NR1STM
            IF (ABS(VVTY(I,K)).LT.EPS12) THEN
              PPLNX(I,K)=0.
              PPLNY(I,K)=1.
            ELSE
              PPLNX(I,K)=1.
              PPLNY(I,K)=-VVTX(I,K)/VVTY(I,K)
            ENDIF
225       CONTINUE
          DO 224 I=1,NR1STM
            PLABS=SQRT(PPLNX(I,K)**2+PPLNY(I,K)**2)
            PPLNX(I,K)=PPLNX(I,K)/PLABS
            PPLNY(I,K)=PPLNY(I,K)/PLABS
224     CONTINUE
C
        DO 227 I=1,NR1STM
        DO 227 J=1,NPPLG
          DO 227 K=NPOINT(1,J),NPOINT(2,J)
            KUP=K+1
            KDN=K
            IF (KUP.GT.NPOINT(2,J)) THEN
              KUP=K
              KDN=K-1
            ENDIF
226         XD=XPOL(I+1,KUP)-XPOL(I+1,KDN)
            YD=YPOL(I+1,KUP)-YPOL(I+1,KDN)
            IF (XD*XD+YD*YD.LT.EPS30) THEN
              IF (KUP.LT.NPOINT(2,J)) THEN
                KUP=KUP+1
              ELSE
                KDN=KDN-1
              ENDIF
              GOTO 226
            ENDIF
            XS=SIGN(1._DP,XD*PPLNX(I,K)+YD*PPLNY(I,K))
            PPLNX(I,K)=PPLNX(I,K)*XS
            PPLNY(I,K)=PPLNY(I,K)*XS
227     CONTINUE
C
C  IDENTIFY DEAD CELLS IN GRID CUTS
C
        IT=1
        DO 229 IR=1,NR1ST-1
        DO 229 J=1,NPPLG-1
          DO 229 IP=NPOINT(2,J),NPOINT(1,J+1)-1
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            NSTGRD(NCELL)=2
229     CONTINUE
C
      ELSEIF (LEVGEO.GT.3) THEN
C
        WRITE (iunout,*) 'ERROR EXIT FROM GRID. NLPOL ',LEVGEO
      ENDIF
C
C  1ST AND 2ND GRID DEFINED
C  SET SURFACE AREA OF NON DEFAULT STANDARD SURFACES
C
C  RADIAL (1ST GRID) SURFACES. OVERWRITE EARLIER VALUES FROM
C  CALL GRID(1)
C
      IF (LEVGEO.EQ.1) THEN
        DO 280 ISTS=1,NSTSI
          IF (INUMP(ISTS,1).NE.0) THEN
            IR=INUMP(ISTS,1)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=PSURF(IRPTE(ISTS,2))-PSURF(IRPTA(ISTS,2))
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              YDF=PSURF(IRPTE(ISTS,2))-PSURF(IRPTA(ISTS,2))
              SAREA(NLJ)=YDF*ZDF*(RSURF(IR)+RMTOR)
            ELSEIF (NLTRT) THEN
            ENDIF
          ENDIF
280     CONTINUE
      ELSEIF (LEVGEO.EQ.2.OR.LEVGEO.EQ.3) THEN
        DO 290 ISTS=1,NSTSI
          IF (INUMP(ISTS,1).NE.0) THEN
            IR=INUMP(ISTS,1)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=BGL(IR,IRPTE(ISTS,2))-BGL(IR,IRPTA(ISTS,2))
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=0.
              DO 291 IP=IRPTA(ISTS,2),IRPTE(ISTS,2)-1
                XS=((XPOL(IR,IP+1)+XPOL(IR,IP))*0.5)+RMTOR
                SAREA(NLJ)=SAREA(NLJ)+(BGL(IR,IP+1)-BGL(IR,IP))*XS
291           CONTINUE
              SAREA(NLJ)=SAREA(NLJ)*TANAL/ALPHA*PI2A
            ENDIF
          ENDIF
290     CONTINUE
      ELSE
C TO BE WRITTEN
      ENDIF
C
C  POLOIDAL (2ND GRID) SURFACES.
C
      IF (LEVGEO.EQ.1) THEN
        DO 285 ISTS=1,NSTSI
          IF (INUMP(ISTS,2).NE.0) THEN
            IP=INUMP(ISTS,2)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=RSURF(IRPTE(ISTS,1))-RSURF(IRPTA(ISTS,1))
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=0.
              DO 286 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(RSURF(IR+1)+RSURF(IR))*0.5+RMTOR
                SAREA(NLJ)=SAREA(NLJ)+(RSURF(IR+1)-RSURF(IR))*XS
286           CONTINUE
              SAREA(NLJ)=SAREA(NLJ)*TANAL/ALPHA*PI2A
            ELSEIF (NLTRT) THEN
              SAREA(NLJ)=0.
              DO 287 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(RSURF(IR+1)+RSURF(IR))*0.5*2.*PIA
                SAREA(NLJ)=SAREA(NLJ)+(RSURF(IR+1)-RSURF(IR))*XS
287           CONTINUE
            ENDIF
          ENDIF
285     CONTINUE
      ELSEIF (LEVGEO.EQ.2.OR.LEVGEO.EQ.3) THEN
        DO 295 ISTS=1,NSTSI
          IF (INUMP(ISTS,2).NE.0) THEN
            IP=INUMP(ISTS,2)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=BGLP(IRPTE(ISTS,1),IP)-BGLP(IRPTA(ISTS,1),IP)
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=0.
              DO 296 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(XPOL(IR+1,IP)+XPOL(IR,IP))*0.5+RMTOR
                SAREA(NLJ)=SAREA(NLJ)+(BGLP(IR+1,IP)-BGLP(IR,IP))*XS
296           CONTINUE
              SAREA(NLJ)=SAREA(NLJ)*TANAL/ALPHA*PI2A
            ELSEIF (NLTRT) THEN
              SAREA(NLJ)=0.
              DO 297 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(XPOL(IR+1,IP)+XPOL(IR,IP))*0.5*2.*PIA
                SAREA(NLJ)=SAREA(NLJ)+(BGLP(IR+1,IP)-BGLP(IR,IP))*XS
297           CONTINUE
            ENDIF
          ENDIF
295     CONTINUE
      ELSE
C TO BE WRITTEN
      ENDIF

!  SET NSTGRD FOR AVERAGING CELLS

      IP = NP2ND
      DO IT = 1, NT3RD
        DO IR = 1, NR1ST
          NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          NSTGRD(NCELL) = 3
        END DO
      END DO

C
      RETURN
C
C   TOROIDAL OR Z-GRID
C
300   CONTINUE
C
C  IF NLSYMT, Z-GRID MUST BE SYMMETRIC
C
      IF (NLTRZ) THEN
C   Z-GRID
        IF (INDGRD(IND).LE.4) THEN
          CALL GRID_1(ZSURF,NT3RD,NTSEP,NTTRA,ZIA,ZGA,ZAA,ZZA,3)
C       ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE Z GRID DATA FROM USER SUPPLIED SUBROUTINE
C TO BE WRITTEN
C         CALL PROUSR (ZSURF,2+4*NPLS+3,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._D
C       ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE Z GRID DATA FROM INTERFACE
C TO BE WRITTEN
C         CALL PROFR (ZSURF,2+4*NPLS+5,1,1,NT3RD)
        ENDIF
C
        DO 310 J=1,NT3RDM
310       ZZONE(J)=(ZSURF(J)+ZSURF(J+1))/2.
C
C     ELSEIF (NLTRA) THEN
C   GRID FOR TOROIDAL APPROXIMATION OF CYLINDER: ALREADY DONE
C
      ENDIF
C
      IF (TRCGRD) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'GRIDPOINTS IN Z DIRECTION'
        CALL LEER(1)
        CALL MASRR1 (' N,  ZSURF ',ZSURF,NT3RD,3)
        CALL LEER(2)
      ENDIF


!  SET NSTGRD FOR AVERAGING CELLS

      IT = NT3RD
      DO IR = 1, NR1ST
        DO IP = 1, NP2ND
          NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          NSTGRD(NCELL) = 3
        END DO
      END DO

!  COPY SWITCHING OFF OF DEAD CELLS FOR TOROIDAL CELL 1 TO ALL OTHERS

      IT1=1
      DO IR=1,NR1ST
        DO IP=1,NP2ND
          NCELL1=IR+((IP-1)+(IT1-1)*NP2T3)*NR1P2
          DO IT=2,NT3RD-1
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            NSTGRD(NCELL)=NSTGRD(NCELL1)
          END DO
        END DO
      END DO
C
      RETURN
C
991   CONTINUE
      WRITE (iunout,*) 'GRID DATA INCONSISTENCY: 3RD GRID.  ZAA > ZIA ?'
      WRITE (iunout,*) 'ZIA,ZAA,NTTRA,ROA= ',ZIA,ZAA,NTTRA,ROA
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'GRID DATA INCONSISTENCY: 2ND GRID.  YAA > YIA ?'
      WRITE (iunout,*) 'YIA,YAA = ',YIA,YAA
      CALL EXIT_OWN(1)
993   CONTINUE
      WRITE (iunout,*) 'GRID DATA INCONSISTENCY: 1ST GRID.  RAA > RIA ?'
      WRITE (iunout,*) 'RIA,RAA = ',RIA,RAA
      WRITE (iunout,*) 'RIA,RAA = ',RIA,RAA
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: input.f
!pb  09.10.06:  save NZADD for higher timesteps
!dr  20.04.06:  fort.10 added as density-model in block 5.
!dr             Also other density model may now refere to test-
!dr             particle tallies, for post processing and iteration
!pb  02.03.06:  NLRAY: switch on raytracing method for stratum 
!pb  20.01.06:  line of sight for cell based spectrum introduced 
!pb  12.01.06:  flag for cell based spectrum added in block 10F  
C    24.09 05:  CALL TO IF2COP NOW LATER, AFTER ALL PLASMA DATA ARE SET
C               SO THAT IF2COP CAN BE USED WITHOUT HAVING TO USE IF1COP
C 
cdr june-05:  spectrum input (10F) extended: SPC_SHIFT,.....
cdr                                SPCPLT_X,SPCPLT_Y,SPCPLT_SAME
cdr           see corresponding changes in CESTIM (ESTIML...)
cdr  28.4.04: nhsts(ispz) introduce, to select species
cdr           for trajectory plot
cdr           default: = 0: "plot trajectory for this species"
cdr           new    : =-1: "do not plot trajectory for this species"
cdr           see modifications in plt2d.f from 28.4.04
cpb sept-05:  specification of filenames for reaction databases in input 
cpb           block 4 added. Lines are inserted between number of reactions
cpb           and reaction cards
cpb           example:   
cpb           CFILE AMJUEL /home/boerner/Database/AMdata/amjuel.tex

      SUBROUTINE INPUT
C
C   READ INPUT DATA AND SET DEFAULT VALUES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CADGEO
      USE CCONA
      USE CGRPTL
      USE CLOGAU
      USE CPL3D
      USE CPLOT
      USE CINIT
      USE COMSIG
      USE CREF
      USE CREFMOD
      USE CPOLYG
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CCOUPL
      USE CGEOM
      USE CSDVI
      USE CSDVI_COP
      USE COMPRT
      USE CPES
      USE COMNNL
      USE COMSOU
      USE CSTEP
      USE COMSPL
      USE CTEXT
      USE CLGIN
      USE COMXS
      USE CSPEI
      USE CTRIG
      USE CTETRA
      USE CESTIM
      USE CUPD
      USE PHOTON

      IMPLICIT NONE

      INCLUDE 'mpif.h'
C
      TYPE TEMPERATURE
        DOUBLE PRECISION          :: TE, TI
        INTEGER                   :: IN, IDION
        TYPE(TEMPERATURE),POINTER :: NEXT
      END TYPE TEMPERATURE
C
      TYPE DENSITY
        DOUBLE PRECISION      :: DI
        INTEGER               :: IN, IDION
        TYPE(DENSITY),POINTER :: NEXT
      END TYPE DENSITY
C
      TYPE VELOCITY
        DOUBLE PRECISION       :: VX, VY, VZ
        INTEGER                :: IZ, IN, IDION
        TYPE(VELOCITY),POINTER :: NEXT
      END TYPE VELOCITY
C
      TYPE VOLUMEP
        DOUBLE PRECISION     :: VOL
        INTEGER              :: IN
        TYPE(VOLUMEP),POINTER :: NEXT
      END TYPE VOLUMEP
C
      TYPE(TEMPERATURE),POINTER :: TEMPLIST, TEMPCUR
      TYPE(DENSITY),POINTER :: DENLIST, DENCUR
      TYPE(VELOCITY),POINTER :: VELLIST, VELCUR
      TYPE(VOLUMEP),POINTER :: VOLLIST, VOLCUR
C
      TYPE SURFACE
        CHARACTER(70)         :: MODNAME
        INTEGER               :: NOSURF
        TYPE(SURFACE),POINTER :: NEXT
      END TYPE SURFACE
C
      TYPE REFMODEL
        CHARACTER(70) :: REFNAME
        INTEGER       :: JLREF,JLSPT
        INTEGER, DIMENSION(:), POINTER :: JSRS,JSRC
        REAL(DP)      :: ZNMLR,EWALLR,EWBINR,FSHEATR
        REAL(DP), DIMENSION(:,:), POINTER :: TRANSPR
        REAL(DP), DIMENSION(:), POINTER ::
     .                                   RCYCFR,RCYCTR,RCPRMR,
     .                                   EXPPLR,EXPELR,EXPILR,
     .                                   RCYCSR,RCYCCR,STPRMR
        TYPE(REFMODEL),POINTER :: NEXT
      END TYPE REFMODEL
C
      TYPE(SURFACE), POINTER :: SURFLIST, SURFCUR, SURFCUR2
      TYPE(REFMODEL), POINTER :: REFLIST, REFCUR

      TYPE REFFILE
        CHARACTER(72) :: RFILE
        TYPE(REFFILE), POINTER :: NEXT
      END TYPE REFFILE

      TYPE(REFFILE), POINTER :: REFFILES, CURFILE

      TYPE(EIRENE_SPECTRUM), POINTER :: ESPEC, SSPEC
C
      REAL(DP) :: AFF(3,3), AFFI(3,3), FP(6)
      REAL(DP) :: RP1, SA, SI, THMAX, SM, SPP, DTIMVO, SAVE, VOLTOT_TAL,
     .          XLREF, YLREF, ZLREF, XLROT, YLROT, ZLROT, XLCOR, SPH,
     .          YLCOR, ZLCOR, VL3, VL4, VL5, ALR, ROTNRM, XSH, RPSDL,
     .          YSH, ALROT, REFNRM, ZSH, VL0, VL1, VL2, RMN, DPP, RMX,
     .          SPCMN, SPCMX,SPC_SHIFT,
     .          SPCPLT_X,SPCPLT_Y,SPCPLT_SAME, SPCVX, SPCVY, SPCVZ, 
     .          VNORM, RCMIN, RCMAX
      REAL(DP) :: tpb1, tpb2, SECOND_OWN
      INTEGER :: IHELP(NLIMPS), IPRSF(12), NUMTAL(12), IADTYP(0:4) 
      INTEGER :: JP, KT, IM, II, IP, IA, NFLGS, NSPZS1, NSPZS2, IPH,
     .           NTLVF, NSRF, NTLS, I1000, NSP, NTL, ICHORI, IRAD,
     .           ILIMPS, ISS, ILA, IB, ILE, INT, NSOPT, IIN, ITEND,
     .           ILLZ, IZ, ISTREAM, ITALI, IAN, IAB, IBEND, IEN, NO,
     .           IGO, IRPTA3, IRPTE2, IRPTE3, ITINI, IH, IDIMP,
     .           JDUMMY, IRPTA1, IRPTA2, IRPTE1, NLJ, I1, I2, I3,
     .           NTIME0, NITER0, IERROR, IREAD, I, ISTS, J, IST,
     .           IPOS2, IPOS0, NM, K, JJ, IPOS1, NRGEN, INUM, NTLSF,
     .           L, INILGJ, INI, ICO, IS, NTLV, ID, IRE,
     .           IN, INELGJ, NPRCSF, MXL, NSPZV1, NSPZV2, NFLGV,
     .           IPRCSF, IR, MT, MP, NDUMM, NUMSEC, NDUMM1, NDUMM2,
     .           NRTAL1, NCOPI, NCOPII, NCOPIE, NFR, NREAC_ADD, IPLN,
     .           NDUMM3, NDUMM4, NRE, ISPSRF, ISPTYP, NSPS, IPTYP, IPSP,
     .           IANF, IEND, IDEFLT_SPUT, IDEFLT_SPEZ, ITLVOUT, NTLVOUT,
     .           ITLSOUT, NTLSOUT, IPLSTI, IPLSV, IFILE, ISRFCLL, 
     .           IDIREC, ISTCHR, JFEXMN, JFEXMX
      INTEGER, SAVE :: NZADD
      INTEGER, EXTERNAL :: IDEZ
      LOGICAL :: LHELP(NLIMPS)
      LOGICAL :: LRPS3D, LRPSCN
      CHARACTER(10) :: CDATE, CTIME
      CHARACTER(80) :: ZEILE, FILE
      CHARACTER(8) :: FILNAM, varname, spcname
      CHARACTER(4) :: H123
      CHARACTER(9) :: REAC
      CHARACTER(50) :: REAC2
      CHARACTER(3) :: CRC
      CHARACTER(60) :: PATH
      CHARACTER(72), ALLOCATABLE, SAVE ::
     .               TXTTLA(:), TXTTLC(:), TXTTLR(:), TXTTLT(:)
      CHARACTER(24), ALLOCATABLE, SAVE ::
     .               TXTSCA(:), TXTSCC(:), TXTSCR(:), TXTSCT(:),
     .               TXTUTA(:), TXTUTC(:), TXTUTR(:), TXTUTT(:)
      CHARACTER(6) :: HANDLE
      CHARACTER(2) :: ELNAME
C
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST ITERATION
C  STEP IN THIS RUN. IITER IS THE ACTUAL ITERATION NUMBER
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST TIMESTEP
C  IN THIS RUN. ITIMV IS THE ACTUAL TIMESTEP NUMBER
C
C  INITIALIZE SOME DATA AND SET DEFAULTS
C
      IREAD=1
C
C  UNIT NUMBER FOR INPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15
      IUNIN=1
C
C  UNIT NUMBER FOR OUTPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15 AND IUNIN
      IUNOUT=6
C
      IF (IITER.GT.1) GOTO 4000
      IF (ITIMV.GT.1) GOTO 4000
      
!pb      TPB1=SECOND_OWN()      
C
      IREAD=0
      IERROR=0
C
      NAINI=0
      NCPVI=0
      NBGVI=0
      MTSURF=0
C
C  SET DEFAULT REACTION MODELS
C
      CALL SETUP_DEFAULT_REACTIONS

C
C  SET DEFAULT SOURCE MODEL
C
      NSTRAI=0
C
C  SET DEFAULT 'ADDITIONAL SURFACE' AND 'STANDARD SURFACE' DATA
C
      NBITS=BIT_SIZE(I)
      CALL SET_DEF_SURF_DATA

      NULLIFY(SURFLIST)
      NULLIFY(REFLIST)
C
C  SET DEFAULT DATA FOR BLOCK 13
C
      DTIMV=1.D30
      TIME0=0.
      NSNVI=0
      NTMSTP=1

C  BY DEFAULT SWITCH OFF MOMENTUM DENSITY TALLIES 
C  FOR USE OF THOSE TALLIES THEY NEED TO BE 
C  SWITCHED ON IN BLOCK 11 EXPLICITELY

C  LV?DEN.. IS AN ALIAS FOR AN ENTRY IN ARRAY LMISTALV
C  THEREFORE .TRUE. MEANS: SWITCHED OFF
      LVXDENA  = .TRUE.
      LVXDENM  = .TRUE.
      LVXDENI  = .TRUE.
      LVXDENPH = .TRUE.
      LVYDENA  = .TRUE.
      LVYDENM  = .TRUE.
      LVYDENI  = .TRUE.
      LVYDENPH = .TRUE.
      LVZDENA  = .TRUE.
      LVZDENM  = .TRUE.
      LVZDENI  = .TRUE.
      LVZDENPH = .TRUE.
C
      CALL LEER(2)
C
99    CONTINUE
C
C  READ TEXT DESCRIBING THE RUN, 100--199
C
100   CONTINUE
C
!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time vo einlesen ',tpb2-tpb1
!pb      tpb1 = tpb2

      CALL DATE_AND_TIME(CDATE,CTIME)
      READ(CDATE(1:4),*) I1
      READ(CDATE(5:6),*) I2
      READ(CDATE(7:8),*) I3
      WRITE (iunout,'(1X,A6,1X,2(I2,1X),I4)') 'DATE: ',I3,I2,I1
      READ(CTIME(1:2),*) I1
      READ(CTIME(3:4),*) I2
      READ(CTIME(5:6),*) I3
      WRITE (iunout,'(1X,A6,1X,3(I2,1X))') 'TIME: ',I1,I2,I3
      CALL LEER(2)
C
      READ (IUNIN,'(A72)') TXTRUN
      WRITE (iunout,'(1X,A72)') TXTRUN
      CALL LEER(1)
109   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1).EQ.'*') THEN
        IF (ZEILE(1:3).NE.'***') WRITE (iunout,'(1X,A72)') ZEILE
        CALL LEER(1)
        GOTO 109
      ELSE
        READ (ZEILE,6666) NMACH,NMODE,NTCPU,NFILE,NITER0,NITER,
     .                    NTIME0,NTIME
      ENDIF
      IITER=MAX0(1,NITER0)
      ITIMV=MAX0(1,NTIME0)

      READ (IUNIN,'(A72)') ZEILE
      IF ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') + INDEX(ZEILE,'T') +
     .     INDEX(ZEILE,'t')) == 0) THEN
        READ (ZEILE,6666) NOPTIM,NOPTM1,NGEOM_USR,NCOUP_INPUT,
     .                    NSMSTRA,NSTORAM,NGSTAL,NRTAL1,NREAC_ADD
        READ (IUNIN,6665) NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
     .                    NLERG,LIDENT,LHABER,LRCHRD
      ELSE
C  THESE DEFAULTS HAVE ALREADY BEEN SET IN FIND-PARAM
!       NOPTIM = 1
!       NOPTM1 = 1
!       NGEOM_USR = 0
!       NCOUP_INPUT = 1
!       NSMSTRA = 1
!       NSTORAM = 9
!       NGSTAL = 0
!       NRTAL1 = 0     ONLY FOR FIND-PARAM, NOT USED ANY FURTHER
!       NREAC_ADD = 0  ONLY FOR FIND-PARAM, NOT USED ANY FURTHER
        READ (ZEILE,6665) NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
     .                    NLERG,LIDENT,LHABER,LRCHRD
      END IF
  
      READ (IUNIN,'(A80)') ZEILE
      IREAD=1
      I1 = INDEX(ZEILE,'CFILE')
      DO WHILE (I1 /= 0)
        I2 = VERIFY(ZEILE(I1+5:),' ') + I1 + 4
        I3 = SCAN(ZEILE(I2+1:),' ')
        HANDLE=REPEAT(' ',6)
        HANDLE(1:I3) = ZEILE(I2:I2+I3-1)
        DO IFILE = 1,NDBNAMES
          IF (INDEX(DBHANDLE(IFILE),HANDLE) /= 0) EXIT
        END DO
        IF (IFILE <= NDBNAMES) THEN
          IANF = I2+I3+VERIFY(ZEILE(I2+I3:),' ')-1
          IEND = IANF+SCAN(ZEILE(IANF+1:),' ')-1
          DBFNAME(IFILE)(1:IEND-IANF+1) = ZEILE(IANF:IEND)
        ELSE
          WRITE (IUNOUT,*) ' WRONG NAME FOR DATABASE ENTERED '
          WRITE (IUNOUT,*) ' DATABASE DEFINITION FOR ',HANDLE,
     .                     ' IGNORED '
        END IF  
        READ (IUNIN,'(A80)') ZEILE       
        I1 = INDEX(ZEILE,'CFILE')
      END DO      

      NFILEN=IDEZ(NFILE,1,5)
      NFILEM=IDEZ(NFILE,2,5)
      NFILEL=IDEZ(NFILE,3,5)
      NFILEK=IDEZ(NFILE,4,5)
      NFILEJ=IDEZ(NFILE,5,5)
      CALL LEER(2)
      CALL MASAGE ('*** 1. DATA FOR OPERATING MODE                   ')
      CALL LEER(1)
      IF (NMACH.EQ.1) THEN
        CALL MASAGE ('       EIRENE RUN ON CRAY                      ')
      ELSEIF (NMACH.EQ.2) THEN
        CALL MASAGE ('       EIRENE RUN ON IBM                       ')
      ELSEIF (NMACH.EQ.3) THEN
        CALL MASAGE ('       EIRENE RUN ON FACOM                     ')
      ELSEIF (NMACH.EQ.4) THEN
        CALL MASAGE ('       EIRENE RUN ON VAX                       ')
      ENDIF
      CALL LEER(1)
      IF (NMODE.NE.0) THEN
        CALL MASAGE ('       EIRENE READS ADDITIONAL DATA FROM FILE  ')
        CALL MASAGE ('       INTERFACING ROUTINE INFCOP IS CALLED    ')
        CALL MASAGE ('       AT ENTRIES IF0COP (GEOMETRY)            ')
        CALL MASAGE ('                  IF1COP (BACKGROUND MEDIUM)   ')
        CALL MASAGE ('                  IF2COP (BOUNDARY CONDITIONS) ')
        IF (NMODE.GT.0) THEN
          CALL MASAGE ('       RETURN DATA TO EXTERNAL CODE:         ')
          CALL MASAGE ('       ENTRIES IF3COP AND IF4COP ARE         ')
          CALL MASAGE ('       CALLED AT THE END OF EACH STRATUM     ')
          CALL MASAGE ('       AND AT THE END OF THE RUN, RESP.      ')
        ELSE
          CALL MASAGE ('       NO RETURN OF DATA TO EXTERNAL CODE:   ')
          CALL MASAGE ('       ENTRIES IF3COP AND IF4COP             ')
          CALL MASAGE ('       ARE NOT CALLED                        ')
        ENDIF
        WRITE (iunout,*) '       NMODE= ',NMODE
      ELSE
        CALL MASAGE ('       EIRENE RUN AS STAND ALONE CODE          ')
        CALL MASAGE ('       INTERFACING ROUTINE INFCOP IS NOT CALLED')
      ENDIF
      CALL LEER(1)
      WRITE (iunout,*) '       EIRENE ASSUMES A TOTAL CPUTIME '
      WRITE (iunout,*) '       OF ',NTCPU,' SECONDS'
      CALL LEER(1)
      IF (NFILEN.EQ.1) THEN
        WRITE (iunout,*) '       EIRENE SAVES OUTPUT DATA '
        WRITE (iunout,*) '       ON FILES FT10 AND FT11 AFTER'
        WRITE (iunout,*) 
     .    '       HAVING COMPUTED THE PARTICLE HISTORIES '
      ELSEIF (NFILEN.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS OUTPUT DATA FROM '
        WRITE (iunout,*) 
     .    '       AN EARLIER RUN FROM FILES FT10 AND FT11 '
        WRITE (iunout,*) '       NO NEW HISTORIES ARE COMPUTED       '
      ELSEIF (NFILEN.EQ.6) THEN
        WRITE (iunout,*) '       EIRENE SAVES OUTPUT DATA '
        WRITE (iunout,*) '       ON FILES FT10 AND FT11 AFTER'
        WRITE (iunout,*) 
     .    '       HAVING COMPUTED THE PARTICLE HISTORIES '
        WRITE (iunout,*) '       FOR THE SUM OVER STRATA TALLIES ONLY '
      ELSEIF (NFILEN.EQ.7) THEN
        WRITE (iunout,*) '       EIRENE READS OUTPUT DATA FROM '
        WRITE (iunout,*) 
     .    '       AN EARLIER RUN FROM FILES FT10 AND FT11 '
        WRITE (iunout,*) '       FOR THE SUM OVER STRATA TALLIES ONLY '
        WRITE (iunout,*) '       NO NEW HISTORIES ARE COMPUTED       '
      ENDIF
      CALL LEER(1)
      IF (NFILEM.EQ.1) THEN
        WRITE (iunout,*) '       EIRENE SAVES GEOMETRICAL DATA '
        WRITE (iunout,*) '       ON FILE FT12'
      ELSEIF (NFILEM.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS GEOMETRICAL DATA FROM '
        WRITE (iunout,*) '       AN EARLIER RUN FROM FILE FT12 '
      ENDIF
      CALL LEER(1)
      IF (NFILEL.EQ.1) THEN
        WRITE (iunout,*) '       EIRENE SAVES PLASMA DATA, A&M DATA'
        WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
        WRITE (iunout,*) 
     .    '       ON FILE FT13 AT END OF RUN, I.E., AFTER'
        WRITE (iunout,*) '       LAST TIMESTEP OR ITERATION '
      ELSEIF (NFILEL.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS PLASMA, A&M DATA      '
        WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA FROM'
        WRITE (iunout,*) '       FILE FT13 '
      ELSEIF (NFILEL.EQ.3) THEN
        WRITE (iunout,*) '       EIRENE READS PLASMA, A&M DATA '
        WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA FROM'
        WRITE (iunout,*) '       FILE FT13  AND '
        WRITE (iunout,*) '       SAVES PLASMA DATA, A&M DATA'
        WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
        WRITE (iunout,*) 
     .    '       ON FILE FT13 AT END OF RUN, I.E., AFTER'
        WRITE (iunout,*) '       LAST TIMESTEP OR ITERATION '
      ELSEIF (NFILEL.EQ.4) THEN
        WRITE (iunout,*) '       EIRENE READS PLASMA AND A&M DATA      '
        WRITE (iunout,*) '       FROM FILE FT13 '
        WRITE (iunout,*) '       SOURCE DISTRIBUTION IS OMITTED '
      ELSEIF (NFILEL.EQ.6) THEN
        WRITE (iunout,*) '       EIRENE SAVES PLASMA DATA, A&M DATA'
        WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
        WRITE (iunout,*) 
     .    '       ON FILE FT13 AT END OF RUN, I.E., AFTER'
        WRITE (iunout,*) '       LAST TIMESTEP OR ITERATION '
        WRITE (iunout,*) '       IN XDR FORMAT              '
      ELSEIF (NFILEL.EQ.7) THEN
        WRITE (iunout,*) '       EIRENE READS PLASMA, A&M DATA      '
        WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA FROM'
        WRITE (iunout,*) '       FILE FT13 '
        WRITE (iunout,*) '       IN XDR FORMAT              '
      ELSEIF (NFILEL.EQ.8) THEN
        WRITE (iunout,*) '       EIRENE READS PLASMA, A&M DATA '
        WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA FROM'
        WRITE (iunout,*) '       FILE FT13  AND '
        WRITE (iunout,*) '       SAVES PLASMA DATA, A&M DATA'
        WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
        WRITE (iunout,*) 
     .    '       ON FILE FT13 AT END OF RUN, I.E., AFTER'
        WRITE (iunout,*) '       LAST TIMESTEP OR ITERATION '
        WRITE (iunout,*) '       IN XDR FORMAT              '
      ELSEIF (NFILEL.EQ.9) THEN
        WRITE (iunout,*) '       EIRENE READS PLASMA AND A&M DATA      '
        WRITE (iunout,*) '       FROM FILE FT13 '
        WRITE (iunout,*) '       SOURCE DISTRIBUTION IS OMITTED '
        WRITE (iunout,*) '       IN XDR FORMAT              '
      ENDIF
      CALL LEER(1)
      IF (NFILEK.EQ.1) THEN
        WRITE (iunout,*) 
     .    '       EIRENE SAVES DATA FOR RECOMMENDED INPUT'
        WRITE (iunout,*) '       MODIFICATIONS ON FILE FT14'
      ELSEIF (NFILEK.EQ.2) THEN
        WRITE (iunout,*) 
     .    '       EIRENE READS DATA FOR RECOMMENDED INPUT'
        WRITE (iunout,*) 
     .    '       MODIFICATIONS FROM FILE FT14, AND CARRIES'
        WRITE (iunout,*) '       THEM OUT IN THIS RUN'
      ELSEIF (NFILEK.EQ.3) THEN
        WRITE (iunout,*) 
     .    '       EIRENE READS OLD DATA FOR RECOMMENDED INPUT'
        WRITE (iunout,*) 
     .    '       MODIFICATIONS FROM FILE FT14, AND CARRIES'
        WRITE (iunout,*) '       THEM OUT IN THIS RUN'
        WRITE (iunout,*) 
     .    '       EIRENE SAVES NEW DATA FOR RECOMMENDED INPUT'
        WRITE (iunout,*) 
     .    '       MODIFICATIONS ON FILE FT14 FOR NEXT RUN'
      ENDIF
      CALL LEER(1)
      IF (NFILEJ.EQ.1) THEN
        WRITE (iunout,*) '       EIRENE SAVES SNAPSHOT POPULATION AT '
        WRITE (iunout,*) '       END OF LAST TIMESTEP ON FILE FT15'
      ELSEIF (NFILEJ.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS SNAPSHOT POPULATION FOR'
        WRITE (iunout,*) '       STRATUM NSTRAI+1 FOR FIRST TIMESTEP'
        WRITE (iunout,*) '       FROM FILE FT15'
      ELSEIF (NFILEJ.EQ.3) THEN
        WRITE (iunout,*) '       EIRENE READS SNAPSHOT POPULATION FOR'
        WRITE (iunout,*) '       STRATUM NSTRAI+1 FOR FIRST TIMESTEP '
        WRITE (iunout,*) '       FROM  FILE FT15 '
        WRITE (iunout,*) '       EIRENE SAVES NEW SNAPSHOT POPULATION'
        WRITE (iunout,*) '       AT END OF LAST TIMESTEP ON FILE FT15'
      ENDIF
      CALL LEER(1)
      IF (NITER.GE.1) THEN
        WRITE (iunout,*) '       EIRENE RUN IN ITERATIVE MODE.      '
        WRITE (iunout,*) '       ITERATIONS: ',IITER,' TO ',NITER
        WRITE (iunout,*) 
     .    '       SUBROUTINE "MODUSR" IS CALLED AFTER EACH'
        WRITE (iunout,*) '       ITERATION '
      ELSE
        WRITE (iunout,*) '       EIRENE RUN IN NONITERATIVE MODE     '
      ENDIF
      CALL LEER(1)
      IF (NTIME.GE.1) THEN
        WRITE (iunout,*) '       EIRENE RUN IN TIME DEP. MODE.       '
        WRITE (iunout,*) '       TIME CYCLES: ',ITIMV,' TO ',NTIME
        WRITE (iunout,*) 
     .    '       SUBROUTINE "TMSUSR" IS CALLED AFTER EACH'
        WRITE (iunout,*) '       TIME CYCLE '
      ELSE
        WRITE (iunout,*) '       EIRENE RUN IN STATIONARY MODE         '
      ENDIF
      CALL LEER(1)
C
C
C  READ DATA FOR STANDARD MESH, 200---299
C
200   CONTINUE
C
      IF (IREAD == 0) READ (IUNIN,*)
      CALL MASAGE ('*** 2. DATA FOR STANDARD MESH                ')
      CALL LEER(1)
C
201   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') THEN
        WRITE (iunout,'(1X,A72)') ZEILE
        CALL LEER(1)
        GOTO 201
      ENDIF
      READ (ZEILE,6666) (INDGRD(J),J=1,3)
C
C INPUT SUB-BLOCK 2A
C
C  RADIAL MESH
210   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') THEN
        GOTO 210
      ENDIF
      READ (ZEILE,6665) NLRAD
      IF (NLRAD) THEN
C
        READ (IUNIN,6665) NLSLB,NLCRC,NLELL,NLTRI,NLPLG,NLFEM,NLTET,
     .                    NLGEN
c slmod begin
        IF (NLTET) THEN 
          READ (IUNIN,6667) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
        ELSE
          READ (IUNIN,6666) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
        ENDIF
c
c        READ (IUNIN,6666) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
c slmod end
        IF (INDGRD(1).LE.5) THEN
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            READ (IUNIN,6664) RIA,RGA,RAA,RRA
            IF (NLELL.OR.NLTRI) THEN
              READ (IUNIN,6664) EP1IN,EP1OT,EP1CH,EXEP1
              READ (IUNIN,6664) ELLIN,ELLOT,ELLCH,EXELL
              IF (NLTRI) THEN
                READ (IUNIN,6664) TRIIN,TRIOT,TRICH,EXTRI
              ENDIF
            ENDIF
          ENDIF
          IF (NLPLG) THEN
            READ(IUNIN,6664) XPCOR,YPCOR,ZPCOR,PLREFL
            READ (IUNIN,6666) (NPOINT(1,K),NPOINT(2,K),K=1,NPPLG)
            DO 212 I=1,NR1ST
              READ (IUNIN,6664) (XPOL(I,J),YPOL(I,J),    J=1,NRPLG)
212         CONTINUE
            IF (PLREFL.GT.0.D0) NR1ST=NR1ST+1
          ENDIF
          IF (NLFEM) THEN
            READ (IUNIN,'(A72)') ZEILE
            IF ((INDEX(ZEILE,'CASE')==0) .AND.
     .          (INDEX(ZEILE,'case') == 0)) THEN
              NTRII=NR1ST
              READ (ZEILE,6664) XPCOR,YPCOR,ZPCOR
              READ (IUNIN,6666) NRKNOT
              READ (IUNIN,6664) (XTRIAN(I),I=1,NRKNOT)
              READ (IUNIN,6664) (YTRIAN(I),I=1,NRKNOT)
              DO 213 J=1,NTRII
                READ (IUNIN,6666) JJ,NECKE (1,J),NECKE (2,J),NECKE(3,J)
                READ (IUNIN,6666)    NCHBAR(1,J),NSEITE(1,J),INMTI(1,J)
                READ (IUNIN,6666)    NCHBAR(2,J),NSEITE(2,J),INMTI(2,J)
                READ (IUNIN,6666)    NCHBAR(3,J),NSEITE(3,J),INMTI(3,J)
213           CONTINUE
            ELSE
              READ (ZEILE(6:),'(A66)') CASENAME
              CASENAME=ADJUSTL(CASENAME)
              I2=INDEX(CASENAME,' ')
              CALL READ_TRIANG (CASENAME(1:I2))
            END IF
          ENDIF
          IF (NLTET) THEN
c slmod begin
            READ (IUNIN,'(A72)') ZEILE
c slmod end
            READ (ZEILE(6:),'(A66)') CASENAME
            CASENAME=ADJUSTL(CASENAME)
            I2=INDEX(CASENAME,' ')
            CALL READ_TETRA (CASENAME(1:I2))
          END IF
          IF (NLGEN) THEN
            NRGEN=NR1ST
          ENDIF
        ELSEIF (INDGRD(1).EQ.6) THEN
C  IS THERE ONE MORE LINE, OR IS NLPOL THE NEXT VARIABLE
          READ (IUNIN,'(A72)') ZEILE
          IPOS1=INDEX(ZEILE,'T')
          IPOS2=INDEX(ZEILE,'F')
          IF (IPOS1.GT.0.OR.IPOS2.GT.0) THEN
            WRITE (iunout,*) 'ONE INPUT LINE MISSING IN BLOCK 2A '
            WRITE (iunout,*) 'AUTOMATIC CORRECTION PERFORMED '
            READ (ZEILE(IPOS1:IPOS1),'(L1)') NLPOL
            GOTO 222
          ENDIF
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            READ (ZEILE,6664) RIA,RGA,RAA
          ELSEIF (NLPLG) THEN
            READ (ZEILE,6664) XPCOR,YPCOR,ZPCOR
          ELSEIF (NLFEM) THEN
            READ (ZEILE,6664) XPCOR,YPCOR,ZPCOR
          ELSEIF (NLTET) THEN
            READ (ZEILE,6664) XPCOR,YPCOR,ZPCOR
          ENDIF
        ENDIF
      ENDIF
C
C  POLOIDAL MESH
C
C INPUT SUB-BLOCK 2B
C
220   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 220
      READ (ZEILE,6665) NLPOL
C
222   READ (IUNIN,6665) NLPLY,NLPLA,NLPLP
      READ (IUNIN,6666) NP2ND,NPSEP,NPPLA,NPPER
      IF (INDGRD(2).LE.5) THEN
        READ (IUNIN,6664) YIA,YGA,YAA,YYA
      ELSEIF (INDGRD(2).EQ.6) THEN
      ENDIF
C
C  TOROIDAL MESH
C
C INPUT SUB-BLOCK 2C
C
230   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 230
      IREAD=1
      READ (ZEILE,6665) NLTOR
      IREAD=0
C
      READ (IUNIN,6665) NLTRZ,NLTRA,NLTRT
      READ (IUNIN,6666) NT3RD,NTSEP,NTTRA,NTPER
      IF (INDGRD(3).LE.5) THEN
        READ (IUNIN,6664) ZIA,ZGA,ZAA,ZZA,ROA
      ELSEIF (INDGRD(3).EQ.6) THEN
      ENDIF
C
C  MESH MULTIPLICATION
C
C INPUT SUB-BLOCK 2D
C
240   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 240
      IREAD=1
      READ (ZEILE,6665) NLMLT
      IREAD=0
C
      IF (NLMLT) THEN
        READ (IUNIN,6666) NBMLT
        READ (IUNIN,6664) (VOLCOR(NM),NM=1,NBMLT)
        READ (IUNIN,'(A72)') ZEILE
      ELSE
        NBMLT=1
        VOLCOR(1)=1.D0
C  FIND START OF NEXT INPUT BLOCK: 2E. SEARCH FOR *, T OR F
241     READ (IUNIN,'(A72)') ZEILE
        IPOS0=INDEX(ZEILE,'*')
        IPOS1=INDEX(ZEILE,'T')
        IPOS2=INDEX(ZEILE,'F')
        IF (IPOS0.EQ.0.AND.IPOS1.EQ.0.AND.IPOS2.EQ.0) GOTO 241
      ENDIF
      IREAD=1
C
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
C
250   IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
C
C INPUT SUB-BLOCK 2E
C
      IF (ZEILE(1:1) .EQ. '*') THEN
        IREAD=0
        GOTO 250
      ENDIF
      READ (ZEILE,6665) NLADD
      IREAD=0
C
      IF (NLADD) THEN
        READ (IUNIN,6666) NRADD
        READ (IUNIN,6664) (VOLADD(NM),NM=1,NRADD)
      ELSE
C  FIND START OF NEXT INPUT BLOCK: 3A
252     READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3) .NE. '***') GOTO 252
        IREAD=1
      ENDIF
C
C
C  READING FOR INPUT BLOCK 2 DONE
C
      NBMLT=MAX0(1,NBMLT)
      NR1ST=MAX0(1,NR1ST)
      NP2ND=MAX0(1,NP2ND)
      IF (.NOT.NLPOL) NP2ND=1
      NT3RD=MAX0(1,NT3RD)
      IF (.NOT.NLTOR) NT3RD=1
C
      IF (NLPLG.AND.INDGRD(1).LE.4) THEN
        IF (NLPOL.AND.NRPLG.NE.NP2ND) GOTO 994
      ENDIF
!pb      IF (NLFEM.AND.INDGRD(1).LT.6) THEN
!pb        GOTO 992
!pb      ENDIF
C
C  READ DATA FOR NON DEFAULT SURFACE MODELS ON STANDARD SURFACES
C  300--349
C
300   CONTINUE
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL MASAGE('*** 3A. DATA FOR NON DEFAULT STANDARD SURFACES')
310   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 310
      IREAD=1
      READ(ZEILE,6666) NSTSI
      IREAD=0
      WRITE (iunout,*) '        NSTSI= ',NSTSI
      CALL LEER(1)
      DO 311 ISTS=1,NSTSI
        NLJ=NLIM+ISTS
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') TXTSFL(NLJ)
        ELSE
          READ (ZEILE,'(A72)') TXTSFL(NLJ)
          IREAD=0
        ENDIF
        READ (IUNIN,6666) JDUMMY,IDIMP,INUMP(ISTS,IDIMP),IRPTA1,
     .                    IRPTE1,IRPTA2,IRPTE2,IRPTA3,IRPTE3

        IF ((IDIMP == 1) .AND. (INUMP(ISTS,IDIMP) > N1ST)) THEN
          WRITE (iunout,*) ' ERROR IN SPECIFICATION OF NON DEFAULT '
          WRITE (iunout,*) ' SURFACE ',ISTS
          WRITE (iunout,*) ' NUMBER OF RADIAL SURFACE > N1ST '
          WRITE (iunout,*) ' CHECK INPUT FILE '
          CALL EXIT_OWN(1)
        ELSEIF ((IDIMP == 2) .AND. (INUMP(ISTS,IDIMP) > N2ND)) THEN
          WRITE (iunout,*) ' ERROR IN SPECIFICATION OF NON DEFAULT '
          WRITE (iunout,*) ' SURFACE ',ISTS
          WRITE (iunout,*) ' NUMBER OF POLOIDAL SURFACE > N2ND '
          WRITE (iunout,*) ' CHECK INPUT FILE '
          CALL EXIT_OWN(1)
        ELSEIF ((IDIMP == 3) .AND. (INUMP(ISTS,IDIMP) > N3RD)) THEN
          WRITE (iunout,*) ' ERROR IN SPECIFICATION OF NON DEFAULT '
          WRITE (iunout,*) ' SURFACE ',ISTS
          WRITE (iunout,*) ' NUMBER OF TOROIDAL SURFACE > N3RD '
          WRITE (iunout,*) ' CHECK INPUT FILE '
          CALL EXIT_OWN(1)
        END IF
C
C  OLD INPUT VERSION BEGIN
        IF (IDIMP.EQ.1.AND.IRPTA1.NE.IRPTE1) THEN
          WRITE (iunout,*) 'WARNING FROM INPUT BLOCK 3A, ISTS= ',ISTS
          WRITE (iunout,*) 'NEW INPUT FOR IRPTA,IRPTE....'
          WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT '
          IRPTA2=IRPTA1
          IRPTE2=IRPTE1
          IRPTA1=INUMP(ISTS,1)
          IRPTE1=INUMP(ISTS,1)
        ELSEIF (IDIMP.EQ.1.AND.IRPTA1.EQ.0.AND.IRPTE1.EQ.0) THEN
          IRPTA1=INUMP(ISTS,1)
          IRPTE1=INUMP(ISTS,1)
        ELSEIF (IDIMP.EQ.2.AND.IRPTA2.EQ.0.AND.IRPTE2.EQ.0) THEN
          IRPTA2=INUMP(ISTS,2)
          IRPTE2=INUMP(ISTS,2)
        ELSEIF (IDIMP.EQ.3.AND.IRPTA3.EQ.0.AND.IRPTE3.EQ.0) THEN
          IRPTA3=INUMP(ISTS,3)
          IRPTE3=INUMP(ISTS,3)
        ENDIF
C  OLD INPUT VERSION DONE
C
C  OVERWRITE DEFAULTS FOR IRPTA, IRPTE ARRAYS
        IF (IDIMP.NE.1) THEN
          IF (IRPTA1.GT.1) IRPTA(ISTS,1)=IRPTA1
          IF (IRPTE1.LE.1.OR.IRPTE1.GT.NR1ST) THEN
            IRPTE(ISTS,1)=MAX(2,NR1ST)
          ELSE
            IRPTE(ISTS,1)=IRPTE1
          ENDIF
          IF ((IDIMP.EQ.2.AND.IRPTE2.GT.NP2ND).OR.
     .        (IDIMP.EQ.3.AND.IRPTE3.GT.NT3RD)) THEN
            WRITE (iunout,*) 'WARNING:'
            WRITE (iunout,*) 'NON-DEFAULT STANDART SURFACE NO. ',ISTS
            WRITE (iunout,*) 'OUT OF RANGE. INVISIBLE!'
          ENDIF
        ENDIF
        IF (IDIMP.NE.2) THEN
          IF (IRPTA2.GT.1) IRPTA(ISTS,2)=IRPTA2
          IF (IRPTE2.LE.1.OR.IRPTE2.GT.NP2ND) THEN
            IRPTE(ISTS,2)=MAX(2,NP2ND)
          ELSE
            IRPTE(ISTS,2)=IRPTE2
          ENDIF
          IF ((IDIMP.EQ.1.AND.IRPTE1.GT.NR1ST).OR.
     .        (IDIMP.EQ.3.AND.IRPTE3.GT.NT3RD)) THEN
            WRITE (iunout,*) 'WARNING:'
            WRITE (iunout,*) 'NON-DEFAULT STANDART SURFACE NO. ',ISTS
            WRITE (iunout,*) 'OUT OF RANGE. INVISIBLE!'
          ENDIF
        ENDIF
        IF (IDIMP.NE.3) THEN
          IF (IRPTA3.GT.1) IRPTA(ISTS,3)=IRPTA3
          IF (IRPTE3.LE.1.OR.IRPTE3.GT.NT3RD) THEN
            IRPTE(ISTS,3)=MAX(2,NT3RD)
          ELSE
            IRPTE(ISTS,3)=IRPTE3
          ENDIF
          IF ((IDIMP.EQ.1.AND.IRPTE1.GT.NR1ST).OR.
     .        (IDIMP.EQ.2.AND.IRPTE2.GT.NP2ND)) THEN
            WRITE (iunout,*) 'WARNING:'
            WRITE (iunout,*) 'NON-DEFAULT STANDART SURFACE NO. ',ISTS
            WRITE (iunout,*) 'OUT OF RANGE. INVISIBLE!'
          ENDIF
        ENDIF
        READ (IUNIN,6666) ILIIN(NLJ),ILSIDE(NLJ),ILSWCH(NLJ),
     .                    ILEQUI(NLJ),ILTOR(NLJ),ILCOL(NLJ),
     .                    ILFIT(NLJ),ILCELL(NLJ),ILBOX(NLJ),
     .                    ILPLG(NLJ)
        IF (ABS(ILCOL(NLJ)).EQ.7) THEN
          WRITE (iunout,*) 'COLOUR FLAG ILCOL CHANGED FOR SURFACE NO. ',
     .                      NLJ
          WRITE (iunout,*) 'COLOUR NO. 7 IS RESERVED FOR "NON-ANALOG'
          WRITE (iunout,*) 
     .      'SURFACES" (SPLITTING, R.R., WEIGHT WINDOWS,..)'
          ILCOL(NLJ)=ILCOL(NLJ)-2
        ENDIF
312     READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1).EQ.'*') THEN
C  SKIP READING LOCAL SURFACE INTERACTION MODEL, USE: DEFAULT
          GOTO 314
        ELSEIF (ILIIN(NLJ).LE.0) THEN
          GOTO 312
        ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
          ALLOCATE(SURFCUR)
          SURFCUR%MODNAME = TRIM(ADJUSTL(ZEILE(9:)))
          SURFCUR%NOSURF = NLJ
          SURFCUR%NEXT => SURFLIST
          SURFLIST => SURFCUR
          IREAD=0
          GOTO 314
        ELSE
          READ (ZEILE,6666) ILREF(NLJ),ILSPT(NLJ),ISRS(1,NLJ),
     .                      ISRC(1,NLJ)
          IREAD=0
          READ (IUNIN,6664) ZNML(NLJ),EWALL(NLJ),EWBIN(NLJ),
     .                      TRANSP(1,1,NLJ),TRANSP(1,2,NLJ),FSHEAT(NLJ)
          READ (IUNIN,6664) RECYCF(1,NLJ),RECYCT(1,NLJ),RECPRM(1,NLJ),
     .                      EXPPL(1,NLJ),EXPEL(1,NLJ),EXPIL(1,NLJ)
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
          IF (ZEILE(1:1).NE.'*') THEN
            READ (ZEILE,6664) RECYCS(1,NLJ),RECYCC(1,NLJ),SPTPRM(1,NLJ)
            IREAD=0
          ELSEIF (ILSPT(NLJ).NE.0) THEN
            WRITE (iunout,*) 'WARNING: SPUTTERING AT NON DEF. SURFACE ',
     .                        ISTS
            WRITE (iunout,*) 
     .        'BUT NO PARAMETERS RECYCS, RECYCC ARE READ '
            WRITE (iunout,*) 'DEFAULT MODEL: "NO SPUTTERING" IS USED. '
            WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
            ILSPT(NLJ)=0
          ENDIF
        ENDIF
C
314     CONTINUE
C
311   CONTINUE
C
C     READ DATA FOR ADDITIONAL SURFACES 350--399
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL MASAGE ('*** 3B. DATA FOR ADDITIONAL SURFACES           ')
350   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1).EQ.'*') THEN
        WRITE (iunout,'(1X,A72)') ZEILE
        GOTO 350
      ENDIF
      IREAD=0
      READ (ZEILE,6666) NLIMI
      WRITE (iunout,*) '        NLIMI= ',NLIMI
      CALL LEER(1)
C
      IF (NLIMI.GT.0) THEN
        DO 353 I=1,NLIMI
          IHELP(I)=IGJUM0(I)
353     CONTINUE
351     READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:3).EQ.'CH0') THEN
          IREAD=0
          CALL DEKEY (ZEILE(4:72),IHELP,1,1,1,NLIMPS)
          GOTO 351
        ENDIF
        DO 352 I=1,NLIMI
          IGJUM0(I)=IHELP(I)
352     CONTINUE
      ENDIF
C
      DO 360 I=1,NLIMI
        IF (IREAD.NE.0) THEN
          READ (ZEILE,'(A72)') TXTSFL(I)
        ELSE
          READ (IUNIN,'(A72)') TXTSFL(I)
        ENDIF
        WRITE (iunout,*) TXTSFL(I)
        IREAD=0
        IF (IGJUM0(I).NE.0) THEN
361       READ (IUNIN,'(A72)') ZEILE
          IREAD=1
          IF (ZEILE(1:1).EQ.'*') GOTO 369
          IF (ZEILE(1:9).EQ.'TRANSFORM') GOTO 368
          GOTO 361
        ENDIF
C
C   GENERAL SURFACE DATA
362     READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3).EQ.'CH1') THEN
          IF (NLIMPB >= NLIMPS) THEN
            CALL DEKEY (ZEILE(4:72),IGJUM1,0,NLIMPS,I,NLIMPS)
          ELSE
            CALL DEKEYB (ZEILE(4:72),IGJUM1,0,NLIMPS,I,NLIMPB,NBITS)
          END IF
          GOTO 362
        ELSEIF (ZEILE(1:3).EQ.'CH2') THEN
          IF (NLIMPB >= NLIMPS) THEN
            CALL DEKEY (ZEILE(4:72),IGJUM2,0,NLIMPS,I,NLIMPS)
          ELSE
            CALL DEKEYB (ZEILE(4:72),IGJUM2,0,NLIMPS,I,NLIMPB,NBITS)
          END IF
          GOTO 362
        ELSE
          READ (ZEILE,6664) RLB(I),SAREA(I),RLWMN(I),RLWMX(I)
          IF (SAREA(I).LE.0.D0) SAREA(I)=666.
          READ (IUNIN,6666) ILIIN(I),ILSIDE(I),ILSWCH(I),
     .                      ILEQUI(I),ILTOR(I),ILCOL(I),
     .                      ILFIT(I),ILCELL(I),ILBOX(I),
     .                      ILPLG(I)
          IF (ABS(ILCOL(I)).EQ.7) THEN
            WRITE (iunout,*) 
     .        'COLOUR FLAG ILCOL CHANGED FOR SURFACE NO. ',I
            WRITE (iunout,*) 'COLOUR NO. 7 IS RESERVED FOR "NON-ANALOG'
            WRITE (iunout,*) 
     .        'SURFACES" (SPLITTING, R.R., WEIGHT WINDOWS,..)'
            ILCOL(I)=ILCOL(I)-2
          ENDIF
C  READ SURFACE COEFFICIENTS
          IF (RLB(I).LT.2.) THEN
            READ (IUNIN,6664) A0LM(I),A1LM(I),A2LM(I),A3LM(I),A4LM(I),
     .                        A5LM(I)
            READ (IUNIN,6664) A6LM(I),A7LM(I),A8LM(I),A9LM(I)
          ELSEIF (RLB(I).LT.3.) THEN
            READ (IUNIN,6664) P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),
     .                        P2(3,I)
          ELSEIF (RLB(I).LT.5.) THEN
            READ (IUNIN,6664) P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),
     .                        P2(3,I)
            READ (IUNIN,6664) P3(1,I),P3(2,I),P3(3,I),P4(1,I),P4(2,I),
     .                        P4(3,I)
          ELSEIF (RLB(I).LT.7.) THEN
            READ (IUNIN,6664) P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),
     .                        P2(3,I)
            READ (IUNIN,6664) P3(1,I),P3(2,I),P3(3,I),P4(1,I),P4(2,I),
     .                        P4(3,I)
            READ (IUNIN,6664) P5(1,I),P5(2,I),P5(3,I),P6(1,I),P6(2,I),
     .                        P6(3,I)
          ENDIF
        ENDIF
C  READ BOUNDARY DATA
        IF (RLB(I).GT.0..AND.RLB(I).LT.2.) THEN
          READ (IUNIN,6664) XLIMS1(1,I),YLIMS1(1,I),ZLIMS1(1,I),
     .                      XLIMS2(1,I),YLIMS2(1,I),ZLIMS2(1,I)
        ELSEIF (RLB(I).LE.0.D0) THEN
          IH=-RLB(I)
          ILIN(I)=IDEZ(IH,1,2)
          ISCN(I)=IDEZ(IH,2,2)
          DO 363 J=1,ILIN(I)
            READ (IUNIN,6664) ALIMS(J,I),XLIMS(J,I),YLIMS(J,I),
     .                        ZLIMS(J,I)
363       CONTINUE
          DO 364 J=1,ISCN(I)
            READ (IUNIN,6664) ALIMS0(J,I),XLIMS1(J,I),YLIMS1(J,I),
     .                        ZLIMS1(J,I),XLIMS2(J,I),YLIMS2(J,I)
            READ (IUNIN,6664) ZLIMS2(J,I),XLIMS3(J,I),YLIMS3(J,I),
     .                        ZLIMS3(J,I)
364       CONTINUE
        ENDIF
C
C READ LOCAL SURFACE INTERACTION MODEL FOR SURFACE NO. I
C
367     READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1).EQ.'*') THEN
C  SKIP READING LOCAL SURFACE INTERACTION MODEL, USE: DEFAULT
          GOTO 369
        ELSEIF (ZEILE(1:9).EQ.'TRANSFORM') THEN
C  TRANSFORM BLOCK
          GOTO 368
        ELSEIF (ILIIN(I).LE.0) THEN
C  TRANSPARENT: SKIP READING SURFACE INTERACTION MODEL
          GOTO 367
        ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
          ALLOCATE(SURFCUR)
          SURFCUR%MODNAME = TRIM(ADJUSTL(ZEILE(9:)))
          SURFCUR%NOSURF = I
          SURFCUR%NEXT => SURFLIST
          SURFLIST => SURFCUR
          IREAD=0
          GOTO 368
        ELSE
C  READ LOCAL REFLECTION MODEL
          READ (ZEILE,6666) ILREF(I),ILSPT(I),ISRS(1,I),ISRC(1,I)
          IREAD=0
          READ (IUNIN,6664) ZNML(I),EWALL(I),EWBIN(I),
     .                      TRANSP(1,1,I),TRANSP(1,2,I),FSHEAT(I)
          READ (IUNIN,6664) RECYCF(1,I),RECYCT(1,I),RECPRM(1,I),
     .                      EXPPL(1,I),EXPEL(1,I),EXPIL(1,I)
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
          IF (ZEILE(1:1).NE.'*'.AND.ZEILE(1:9).NE.'TRANSFORM') THEN
            READ (ZEILE,6664) RECYCS(1,I),RECYCC(1,I),SPTPRM(1,I)
            IREAD=0
          ELSEIF (ILSPT(NLJ).NE.0) THEN
            WRITE (iunout,*) 'WARNING: SPUTTERING FOR ADD. SURFACE ',NLJ
            WRITE (iunout,*) 
     .        'BUT NO PARAMETERS RECYCS, RECYCC ARE READ '
            WRITE (iunout,*) 'DEFAULT MODEL: "NO SPUTTERING" IS USED. '
            WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
            ILSPT(NLJ)=0
          ENDIF
        ENDIF
C
368     IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:9).NE.'TRANSFORM') GOTO 369
C
C  TRANSFORM FROM CONVENIENT CO-ORDINATE SYSTEM TO EIRENE CO-ORDINATE
C  SYSTEM. THIS IS POSSIBLE FOR ALL SURFACES READ BY NOW, I.E. FROM
C  ITINI=1 TO ITEND=I
        IREAD=0
        READ (IUNIN,6666) ITINI,ITEND
        IF (ITINI.LT.1) ITINI=I
        IF (ITEND.GT.I) ITEND=I
        READ (IUNIN,6664) XLCOR,YLCOR,ZLCOR
        READ (IUNIN,6664) XLREF,YLREF,ZLREF
        READ (IUNIN,6664) XLROT,YLROT,ZLROT,ALROT
C  SHIFT
        XSH=XLCOR
        IF (XSH.NE.0.D0) CALL XSHADD(XSH,ITINI,ITEND)
        YSH=YLCOR
        IF (YSH.NE.0.D0) CALL YSHADD(YSH,ITINI,ITEND)
        ZSH=ZLCOR
        IF (ZSH.NE.0.D0) CALL ZSHADD(ZSH,ITINI,ITEND)
C  REFLECTION
        REFNRM=XLREF*XLREF+YLREF*YLREF+ZLREF*ZLREF
        IF (REFNRM.GT.EPS10) THEN
          CALL SETREF(AFF,AFFI,1,XLREF,YLREF,ZLREF)
          CALL ROTADD(AFF,AFFI,ITINI,ITEND)
        ENDIF
C  ROTATION
        ROTNRM=XLROT*XLROT+YLROT*YLROT+ZLROT*ZLROT
        ALR=ALROT
        IF (ALR.NE.0..AND.ROTNRM.GT.EPS10) THEN
          CALL SETROT(AFF,AFFI,1,XLROT,YLROT,ZLROT,ALROT)
          CALL ROTADD(AFF,AFFI,ITINI,ITEND)
        ENDIF
        GOTO 368
C
369     CONTINUE
C
360   CONTINUE
C
C  READ DATA FOR SPECIES SPECIFICATION AND ATOMIC PHYSICS MODULE
C  400--499
C
400   CONTINUE
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      IREAD=0
      CALL MASAGE ('*** 4. DATA FOR SPECIES SPECIFICATION AND   ')
      CALL MASAGE ('       ATOMIC PHYSICS MODULE                ')
      CALL LEER(1)
C
      READ (IUNIN,*)
      WRITE (iunout,*) 
     .  '       ATOMIC REACTION CARDS, NREACI DATA FIELDS'
      READ (IUNIN,*) NREACI
      WRITE (iunout,*) '       NREACI= ',NREACI
      CALL LEER(1)
      IF (NPHOTI > 0) CALL PH_INIT(0)
C
411   READ (IUNIN,'(A80)') ZEILE
      IF (ZEILE(1:1).NE.'*') THEN
C
C  READ ONE REACTION FROM FILE "FILNAM"
C
        IF (INDEX(ZEILE,'PHOTON') == 0) THEN
          READ (ZEILE,66661)
     .          IR,FILNAM,H123,REAC,CRC,MP,MT,DPP,RMN,RMX
          REAC2 = REAC
        ELSE
          ZEILE=ADJUSTL(ZEILE)
          I1=SCAN(ZEILE,' ')
          READ (ZEILE(1:I1-1),*) IR
!  READ FILNAM
          IANF=I1+VERIFY(ZEILE(I1+1:),' ')
          IEND=IANF+SCAN(ZEILE(IANF+1:),' ')-1
          IF (IEND-IANF+1 > 8) THEN
            WRITE (iunout,*) ' FILENAME FOR REACTION ',IR,' TOO LONG '
            CALL EXIT_OWN(1)
          END IF
          FILNAM(1:IEND-IANF+1)=ZEILE(IANF:IEND)
!  READ H123
          IANF=IEND+VERIFY(ZEILE(IEND+1:),' ')
          IEND=IANF+SCAN(ZEILE(IANF+1:),' ')-1
          IF (IEND-IANF+1 > 4) THEN
            WRITE (iunout,*) 
     .        ' STRING H123 FOR REACTION ',IR,' TOO LONG '
            CALL EXIT_OWN(1)
          END IF
          H123='    '
          H123(1:IEND-IANF+1)=ZEILE(IANF:IEND)
!  READ REAC
          IANF=IEND+VERIFY(ZEILE(IEND+1:),' ')
          IEND=IANF+SCAN(ZEILE(IANF+1:),' ')-1
          IF (IEND-IANF+1 > 50) THEN
            WRITE (iunout,*) 
     .        ' REACTIONSTRING FOR REACTION ',IR,' TOO LONG '
            CALL EXIT_OWN(1)
          END IF
          REAC2=REPEAT(' ',50)
          REAC2(1:IEND-IANF+1)=ZEILE(IANF:IEND)
!  READ CRC
          IANF=IEND+VERIFY(ZEILE(IEND+1:),' ')
          IEND=IANF+SCAN(ZEILE(IANF+1:),' ')-1
          IF (IEND-IANF+1 > 3) THEN
            WRITE (iunout,*) ' CRC-STRING FOR REACTION ',IR,' TOO LONG '
            CALL EXIT_OWN(1)
          END IF
          CRC(1:IEND-IANF+1)=ZEILE(IANF:IEND)

          READ (ZEILE(IEND+1:),*) MP,MT,DPP,RMN,RMX
        END IF

        IF (INDEX(ZEILE,'ADAS') .NE. 0) THEN
          READ (IUNIN,'(4X,A2,1X,I3)') ELNAME,IZ
          CALL LOWERCASE(ELNAME)
        ELSE
          ELNAME = ' '
          IZ = 0
        END IF

        IF (IR.GT.NREACI) THEN
            CALL MASPRM('NREACI',6,NREACI,'IR',2,IR,IERROR)
        ENDIF
        MASSP(IR)=MP
        MASST(IR)=MT
        DELPOT(IR)=DPP
C  IFEXMN,IFEXMX,FPARM:
C  ASYMPTOTICS FOR CROSS SECTIONS (SECOND INDEX=1) OR RATES (SECOND
C  INDEX=2) , OVERWRITES ASYMPTOTICS IN DATA FILES, IF THERE ARE SUCH
csw added branch
        FP = 0._DP
        IF (INDEX(H123,'P.').eq.0) then
          IF (INDEX(H123,'H.1 ').NE.0) J=1
          IF (INDEX(H123,'H.1 ').EQ.0) J=2
          RCMIN = -20.
          RCMAX =  20.
          JFEXMN = 0
          JFEXMX = 0
          IF (RMN.GT.0.D0) THEN
            READ (IUNIN,66664) JFEXMN,(FP(I),I=1,3)
            RCMIN=LOG(RMN)
          ENDIF
          IF (RMX.GT.0.D0) THEN
            READ (IUNIN,66664) JFEXMX,(FP(I),I=4,6)
            RCMAX=LOG(RMX)
          ENDIF
        else
          RCMIN=-20.
          RCMAX= 20.
          JFEXMN=0
          JFEXMX=0
        endif
csw end branch
C
        CALL SLREAC (IR,FILNAM,H123,REAC2,CRC,
     .               RCMIN,RCMAX,FP,JFEXMN,JFEXMX,ELNAME,IZ)
        GOTO 411
      ENDIF
C
420   WRITE (iunout,*) 
     .  '*4A.   NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES'
      READ (IUNIN,*) NATMI
      WRITE (iunout,*) '       NATMI= ',NATMI
      CALL LEER(1)
C
      NSPH=NPHOTI
      DO 421 IATM=1,NATMI
        ISPZ=NSPH+IATM
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSA(IATM),NCHARA(IATM),
     .                       NDUMM1,NDUMM2,
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCA(IATM),NFOLA(IATM),NGENA(IATM),
     .                       NHSTS(ISPZ)
C  DEFAULTS FOR ATOMIC SPECIES:
        NPRT(ISPZ)=1
        DO 422 K=1,NRCA(IATM)
          IF (NUMSEC .NE. 3) THEN
            READ (IUNIN,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                        ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                        ISCDEA(IATM,K),IESTMA(IATM,K),
     .                        IBGKA(IATM,K)
          ELSE
            READ (IUNIN,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                        ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                        ISCD3A(IATM,K),
     .                        ISCDEA(IATM,K),IESTMA(IATM,K),
     .                        IBGKA(IATM,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*) 
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1A = ',ISCD1A(IATM,K)
            WRITE (iunout,*) ' ISCD2A = ',ISCD2A(IATM,K)
            WRITE (iunout,*) ' ISCD3A = ',ISCD3A(IATM,K)
            WRITE (iunout,*) ' ISCDEA = ',ISCDEA(IATM,K)
            WRITE (iunout,*) ' IESTMA = ',IESTMA(IATM,K)
            WRITE (iunout,*) ' IBGKA  = ',IBGKA(IATM,K)
          END IF
          READ (IUNIN,6664) EELECA(IATM,K),EBULKA(IATM,K),
     .                      ESCD1A(IATM,K),ESCD2A(IATM,K),
     .                      FREACA(IATM,K),FLDLMA(IATM,K)
422     CONTINUE
421   CONTINUE
C
C  READ NEUTRAL MOLECULES SPECIES CARDS
C
      READ (IUNIN,*)
      WRITE (iunout,*) 
     .  '*4B.   NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES'
430   READ (IUNIN,*) NMOLI
      WRITE (iunout,*) '       NMOLI= ',NMOLI
      CALL LEER(1)
      NSPA=NSPH+NATMI
      DO 431 IMOL=1,NMOLI
        ISPZ=NSPA+IMOL
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSM(IMOL),NCHARM(IMOL),
     .                       NPRT(ISPZ),NDUMM,
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCM(IMOL),NFOLM(IMOL),NGENM(IMOL),
     .                       NHSTS(ISPZ)
        IF (ISRT(ISPZ,1).LT.0) THEN
          WRITE (iunout,*) 'INPUT ERROR IN BLOCK 4B '
          WRITE (iunout,*) 'MOLECULAR SPECIES ',IMOL,':'
          WRITE (iunout,*) 'ISRT LT 0 OPTION IS NOT AVAILABLE ANYMORE'
          WRITE (iunout,*) 'PROBABLY YOU MEAN: ISRT= ',NMOLI+1
          IERROR=IERROR+1
        ENDIF
        DO 432 K=1,NRCM(IMOL)
          IF (NUMSEC .NE. 3) THEN
            READ (IUNIN,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                        ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                        ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                        IBGKM(IMOL,K)
          ELSE
            READ (IUNIN,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                        ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                        ISCD3M(IMOL,K),
     .                        ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                        IBGKM(IMOL,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*) 
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1M = ',ISCD1M(IMOL,K)
            WRITE (iunout,*) ' ISCD2M = ',ISCD2M(IMOL,K)
            WRITE (iunout,*) ' ISCD3M = ',ISCD3M(IMOL,K)
            WRITE (iunout,*) ' ISCDEM = ',ISCDEM(IMOL,K)
            WRITE (iunout,*) ' IESTMM = ',IESTMM(IMOL,K)
            WRITE (iunout,*) ' IBGKM  = ',IBGKM(IMOL,K)
          END IF
          READ (IUNIN,6664) EELECM(IMOL,K),EBULKM(IMOL,K),
     .                      ESCD1M(IMOL,K),ESCD2M(IMOL,K),FREACM(IMOL,K)
432     CONTINUE
431   CONTINUE
C
C  READ TEST PARTICLE IONS SPECIES CARDS
C
      READ (IUNIN,*)
      WRITE (iunout,*) '*4C.   TEST IONS SPECIES CARDS, NIONI SPECIES '
440   READ (IUNIN,*) NIONI
      WRITE (iunout,*) '       NIONI= ',NIONI
      CALL LEER(1)
      NSPAM=NSPH+NATMI+NMOLI
      DO 441 IION=1,NIONI
        ISPZ=NSPAM+IION
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSI(IION),NCHARI(IION),
     .                       NPRT(ISPZ),NCHRGI(IION),
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCI(IION),NFOLI(IION),NGENI(IION),
     .                       NHSTS(ISPZ)
        DO 442 K=1,NRCI(IION)
          IF (NUMSEC .NE. 3) THEN
            READ (IUNIN,6666) IREACI(IION,K),IBULKI(IION,K),
     .                        ISCD1I(IION,K),ISCD2I(IION,K),
     .                        ISCDEI(IION,K),IESTMI(IION,K),
     .                        IBGKI(IION,K)
          ELSE
            READ (IUNIN,6666) IREACI(IION,K),IBULKI(IION,K),
     .                        ISCD1I(IION,K),ISCD2I(IION,K),
     .                        ISCD3I(IION,K),
     .                        ISCDEI(IION,K),IESTMI(IION,K),
     .                        IBGKI(IION,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*) ' THREE SECONDARY GROUPS FOR SPECIES ',
     .                    TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1I = ',ISCD1I(IION,K)
            WRITE (iunout,*) ' ISCD2I = ',ISCD2I(IION,K)
            WRITE (iunout,*) ' ISCD3I = ',ISCD3I(IION,K)
            WRITE (iunout,*) ' ISCDEI = ',ISCDEI(IION,K)
            WRITE (iunout,*) ' IESTMI = ',IESTMI(IION,K)
            WRITE (iunout,*) ' IBGKI  = ',IBGKI(IION,K)
          END IF
          READ (IUNIN,6664) EELECI(IION,K),EBULKI(IION,K),
     .                      ESCD1I(IION,K),ESCD2I(IION,K),FREACI(IION,K)
442     CONTINUE
441   CONTINUE
C
      READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3) == '***') THEN
        IREAD=1
        NPHOTI=0
        GOTO 500
      END IF
      IREAD=0
      WRITE (iunout,*) 
     .  '*4D.   NEUTRAL PHOTONS SPECIES CARDS, NPHOTI SPECIES'
      READ (IUNIN,*) NPHOTI
      WRITE (iunout,*) '       NPHOTI= ',NPHOTI
      CALL LEER(1)
C
      DO 451 IPHOT=1,NPHOTI
        ISPZ=IPHOT
        READ (IUNIN,66666) I,TEXTS(ISPZ),NDUMM1,NDUMM2,
     .                       NDUMM3,NDUMM4,
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCPH(IPHOT),NFOLPH(IPHOT),NGENPH(IPHOT),
     .                       NHSTS(ISPZ)
C  DEFAULTS FOR PHOTONIC SPECIES:
        NPRT(ISPZ)=1
        DO 452 K=1,NRCPH(IPHOT)
          IF (NUMSEC .NE. 3) THEN
            READ (IUNIN,6666) IREACPH(IPHOT,K),IBULKPH(IPHOT,K),
     .                        ISCD1PH(IPHOT,K),ISCD2PH(IPHOT,K),
     .                        ISCDEPH(IPHOT,K),IESTMPH(IPHOT,K),
     .                        IBGKPH(IPHOT,K)
          ELSE
            READ (IUNIN,6666) IREACPH(IPHOT,K),IBULKPH(IPHOT,K),
     .                        ISCD1PH(IPHOT,K),ISCD2PH(IPHOT,K),
     .                        ISCD3PH(IPHOT,K),
     .                        ISCDEPH(IPHOT,K),IESTMPH(IPHOT,K),
     .                        IBGKPH(IPHOT,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*) 
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1PH = ',ISCD1PH(IPHOT,K)
            WRITE (iunout,*) ' ISCD2PH = ',ISCD2PH(IPHOT,K)
            WRITE (iunout,*) ' ISCD3PH = ',ISCD3PH(IPHOT,K)
            WRITE (iunout,*) ' ISCDEPH = ',ISCDEPH(IPHOT,K)
            WRITE (iunout,*) ' IESTMPH = ',IESTMPH(IPHOT,K)
            WRITE (iunout,*) ' IBGKPH  = ',IBGKPH(IPHOT,K)
          END IF
          READ (IUNIN,6664) EELECPH(IPHOT,K),EBULKPH(IPHOT,K),
     .                      ESCD1PH(IPHOT,K),ESCD2PH(IPHOT,K),
     .                      FREACPH(IPHOT,K),FLDLMPH(IPHOT,K)
452     CONTINUE
451   CONTINUE
C
C  READ DATA FOR PLASMA-BACKGROUND , 500--599
C
500   CONTINUE
C
      IF (IREAD == 0)  READ (IUNIN,'(A72)') ZEILE
C     WRITE (iunout,*) ZEILE
      CALL MASAGE ('*** 5. DATA FOR PLASMA BACKGROUND            ')
      CALL LEER(1)
C
C  READ BULK IONS SPECIES CARDS
C
      READ (IUNIN,'(A72)') ZEILE
C     WRITE (iunout,*) ZEILE
      WRITE (iunout,*) '*5A.   BULK ION SPECIES CARDS, NPLSI SPECIES '
510   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 510
      READ(ZEILE,6666) NPLSI
      WRITE (iunout,*) '       NPLSI= ',NPLSI
      CALL LEER(1)
      NSPAMI=NSPAM+NIONI
      NSPTOT=NSPAMI+NPLSI
      DO 511 IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSP(IPLS),NCHARP(IPLS),
     .                       NPRT(ISPZ),NCHRGP(IPLS),
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCP(IPLS),NDUMM1,NDUMM2,
     .                       NHSTS(ISPZ),NDUMM4,
     .                       CDENMODEL(IPLS),NRE
        CALL UPPERCASE (CDENMODEL(IPLS))
        DO 512 K=1,NRCP(IPLS)
          IF (NUMSEC .NE. 3) THEN
            READ (IUNIN,6666) IREACP(IPLS,K),IBULKP(IPLS,K),
     .                        ISCD1P(IPLS,K),ISCD2P(IPLS,K),
     .                        ISCDEP(IPLS,K)
          ELSE
            READ (IUNIN,6666) IREACP(IPLS,K),IBULKP(IPLS,K),
     .                        ISCD1P(IPLS,K),ISCD2P(IPLS,K),
     .                        ISCD3P(IPLS,K),ISCDEP(IPLS,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*) 
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1P = ',ISCD1P(IPLS,K)
            WRITE (iunout,*) ' ISCD2P = ',ISCD2P(IPLS,K)
            WRITE (iunout,*) ' ISCD3P = ',ISCD3P(IPLS,K)
            WRITE (iunout,*) ' ISCDEP = ',ISCDEP(IPLS,K)
          END IF
          READ (IUNIN,6664) EELECP(IPLS,K),EBULKP(IPLS,K),
     .                      ESCD1P(IPLS,K),ESCD2P(IPLS,K),FREACP(IPLS,K)
512     CONTINUE
        IF (LEN_TRIM(CDENMODEL(IPLS)) > 0) THEN
          ALLOCATE (TDMPAR(IPLS)%TDM)
          TDMPAR(IPLS)%TDM%NRE=MAX(NRE,1)
          ALLOCATE (TDMPAR(IPLS)%TDM%ISP(TDMPAR(IPLS)%TDM%NRE))
          ALLOCATE (TDMPAR(IPLS)%TDM%ITP(TDMPAR(IPLS)%TDM%NRE))
          ALLOCATE (TDMPAR(IPLS)%TDM%ISTR(TDMPAR(IPLS)%TDM%NRE))
          ALLOCATE (TDMPAR(IPLS)%TDM%FNAME(TDMPAR(IPLS)%TDM%NRE))
          ALLOCATE (TDMPAR(IPLS)%TDM%H2(TDMPAR(IPLS)%TDM%NRE))
          ALLOCATE (TDMPAR(IPLS)%TDM%REACTION(TDMPAR(IPLS)%TDM%NRE))
          ALLOCATE (TDMPAR(IPLS)%TDM%CR(TDMPAR(IPLS)%TDM%NRE))
          SELECT CASE (CDENMODEL(IPLS))
          CASE ('FORT.13   ')
            READ (IUNIN,6666) TDMPAR(IPLS)%TDM%ISP(1) 
                              TDMPAR(IPLS)%TDM%ITP(1)=4
          CASE ('FORT.10   ')
            READ (IUNIN,6666) TDMPAR(IPLS)%TDM%ISP(1), 
     .                        TDMPAR(IPLS)%TDM%ITP(1),
     .                        TDMPAR(IPLS)%TDM%ISTR(1)
          CASE ('SAHA      ')
!PB   TO BE WRITTEN
          CASE ('BOLTZMANN ')
            READ (IUNIN,'(3I6,6x,2E12.4)')
     .           TDMPAR(IPLS)%TDM%ISP(1),
     .           TDMPAR(IPLS)%TDM%ITP(1),
     .           TDMPAR(IPLS)%TDM%ISTR(1),
     .           TDMPAR(IPLS)%TDM%G_BOLTZ, 
     .           TDMPAR(IPLS)%TDM%DELTAE
          CASE ('CORONA    ')
            READ (IUNIN,'(3I6,1X,A6,1X,A4,A9,A3,E12.4)')
     .           TDMPAR(IPLS)%TDM%ISP(1),
     .           TDMPAR(IPLS)%TDM%ITP(1), 
     .           TDMPAR(IPLS)%TDM%ISTR(1),
     .           TDMPAR(IPLS)%TDM%FNAME(1),
     .           TDMPAR(IPLS)%TDM%H2(1),
     .           TDMPAR(IPLS)%TDM%REACTION(1),
     .           TDMPAR(IPLS)%TDM%CR(1),
     .           TDMPAR(IPLS)%TDM%A_CORONA
            IF (INDEX(TDMPAR(IPLS)%TDM%H2(1),'H.2') == 0) THEN
              WRITE (iunout,*) 
     .          ' WRONG REACTION SPECIFIED FOR CORONA MODEL '
              WRITE (iunout,*) ' ONLY H.2 REACTIONS ARE PERMITTED '
              WRITE (iunout,*) ' IPLS = ',IPLS
              CALL EXIT_OWN(1)
            END IF
          CASE ('COLRAD    ')
            DO I=1, TDMPAR(IPLS)%TDM%NRE
              READ (IUNIN,'(3I6,1X,A6,1X,A4,A9,A3)')
     .             TDMPAR(IPLS)%TDM%ISP(I),
     .             TDMPAR(IPLS)%TDM%ITP(I),
     .             TDMPAR(IPLS)%TDM%ISTR(I),
     .             TDMPAR(IPLS)%TDM%FNAME(I),
     .             TDMPAR(IPLS)%TDM%H2(I),
     .             TDMPAR(IPLS)%TDM%REACTION(I),
     .             TDMPAR(IPLS)%TDM%CR(I)
              IF (INDEX(TDMPAR(IPLS)%TDM%H2(1),'H.11') == 0. AND.
     .            INDEX(TDMPAR(IPLS)%TDM%H2(1),'H.12') == 0.) THEN
                WRITE (iunout,*)
     .            ' WRONG REACTION SPECIFIED FOR COLRAD MODEL '
                WRITE (iunout,*) ' ONLY H.11 OR H.12 REACTIONS '
                WRITE (IUNOUT,*) ' ARE PERMITTED '
                WRITE (iunout,*) ' IPLS = ',IPLS
                CALL EXIT_OWN(1)
              END IF
            END DO
          CASE DEFAULT
!PB  NOTHING TO BE DONE
          END SELECT
        END IF
511   CONTINUE
C
      DO I=1,NSPZ
        CALL UPPERCASE (TEXTS(I))
      END DO
C
520   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') THEN
C       WRITE (iunout,*) ZEILE
        GOTO 520
      ENDIF
      CALL MASAGE ('*5B.   PLASMA BACKGROUND DATA                ')
      READ (ZEILE,6666) (INDPRO(J),J=1,12)
C  Te profile
      IF (INDPRO(1).LE.5.AND.NPLSI.GT.0)
     .  READ (IUNIN,6664) TE0,TE1,TE2,TE3,TE4,TE5
C  Ti profile(s)
      NPLSTI = 1
      IF ((INDPRO(2) < 0) .OR. (INDPRO(2) > 9)) NPLSTI=NPLS
      NLMLTI=NPLSTI > 1
      MPLSTI=1
      IF (NLMLTI) MPLSTI = (/ (I,I=1,NPLS) /)
      INDPRO(2)=IABS(INDPRO(2))
      IF (INDPRO(2) > 9) INDPRO(2) = MOD(INDPRO(2),10)
      IF (INDPRO(2).LE.5.AND.NPLSI.GT.0) THEN
        IF (NLMLTI) THEN
          READ (IUNIN,6664) (TI0(I),TI1(I),TI2(I),TI3(I),TI4(I),TI5(I),
     .                       I=1,NPLSI)
        ELSE
C  ONLY ONE COMMON ION TEMPERATURE FOR ALL SPECIES
          READ (IUNIN,6664)  TI0(1),TI1(1),TI2(1),TI3(1),TI4(1),TI5(1)
          DO 530 IPLS=2,NPLSTI
            TI0(IPLS)=TI0(1)
            TI1(IPLS)=TI1(1)
            TI2(IPLS)=TI2(1)
            TI3(IPLS)=TI3(1)
            TI4(IPLS)=TI4(1)
            TI5(IPLS)=TI5(1)
530       CONTINUE
        ENDIF
      ENDIF
c  di profiles
      IF (INDPRO(3).LE.5)
     .  READ (IUNIN,6664) (DI0(I),DI1(I),DI2(I),DI3(I),DI4(I),DI5(I),
     .                     I=1,NPLSI)
c  vi profile(s)
      NLMACH=INDPRO(4).LT.0
      NPLSV = NPLS
      IF (IABS(INDPRO(4)) > 9) NPLSV = 1

      INDPRO(4)=IABS(INDPRO(4))
      IF (INDPRO(4) > 9) INDPRO(4) = MOD(INDPRO(4),10)
      NLMLV = NPLSV > 1
      IF (.not.NLMLV) THEN
        MPLSV = 1
      ELSE
        MPLSV = (/ (I,I=1,NPLS) /)
      ENDIF
      IF (INDPRO(4).LE.5) THEN
        IF (.not. NLMLV) THEN
          READ (IUNIN,6664) VX0(1),VX1(1),VX2(1),VX3(1),VX4(1),VX5(1)
          READ (IUNIN,6664) VY0(1),VY1(1),VY2(1),VY3(1),VY4(1),VY5(1)
          READ (IUNIN,6664) VZ0(1),VZ1(1),VZ2(1),VZ3(1),VZ4(1),VZ5(1)
        ELSE
          READ (IUNIN,6664) (VX0(I),VX1(I),VX2(I),VX3(I),VX4(I),VX5(I),
     .                       I=1,NPLSI)
          READ (IUNIN,6664) (VY0(I),VY1(I),VY2(I),VY3(I),VY4(I),VY5(I),
     .                       I=1,NPLSI)
          READ (IUNIN,6664) (VZ0(I),VZ1(I),VZ2(I),VZ3(I),VZ4(I),VZ5(I),
     .                       I=1,NPLSI)
        END IF
      ENDIF
c  pitch -profile
      IF (INDPRO(5).LE.5)
     .  READ (IUNIN,6664) B0,B1,B2,B3,B4,B5
c  cell volume -profile
      IF (INDPRO(12).LE.5) THEN
        READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:3) .EQ. '***') THEN
          WRITE (iunout,*) 'ONE INPUT LINE MISSING IN BLOCK 5 '
          WRITE (iunout,*) 'AUTOMATIC CORRECTION PERFORMED '
          VL0=0
        ELSE
          READ (ZEILE,6664) VL0,VL1,VL2,VL3,VL4,VL5
          IREAD=0
        ENDIF
      ENDIF
C
C  READ  DATA FOR REFLECTION MODEL  600--699
C
600   CONTINUE
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      IREAD=0
      NULLIFY(REFFILES)
      CALL MASAGE ('*** 6. GENERAL DATA FOR REFLECTION MODEL    ')
      CALL LEER(1)
C
610   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 610
      READ (ZEILE,6665) NLTRIM
      IREAD=0
      IF (NLTRIM) THEN
        READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (INDEX(ZEILE,'PATH')+INDEX(ZEILE,'path').EQ.0) THEN
C  NO PATH SPECIFIED FOR REFLECTION DATA BASE
          WRITE (iunout,*) 
     .      ' NO PATH SPECIFIED FOR REFLECTION DATA BASE '
          WRITE (iunout,*) ' OLD VERSION USED '
          LTRMOL=.TRUE.
C  SKIP LINES CONTAINING SPECIFICATIONS FOR DATA BASES AND
C  CONTINUE READING
615       IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').EQ.0) THEN
            READ (ZEILE,6664) (DATD(IATM),IATM=1,NATMI)
            IREAD=0
            GOTO 620
          ELSE
            READ (IUNIN,'(A72)') ZEILE
            IREAD=1
            GOTO 615
          ENDIF
        ELSE
C  PATH SPECIFICATION FOR DATA BASE FOUND
          LTRMOL=.FALSE.
          READ (ZEILE(7:),'(A60)') PATH
          WRITE (iunout,*) ' PATH = ',PATH
          PATH=ADJUSTL(PATH)
          I2=INDEX(PATH,' ')
          IF ((I2 == 0) .AND. (ZEILE(68:68) /= ' ')) THEN
            WRITE (iunout,*) ' PATH FOR TRIM REFLECTION DATABASE IS',
     .                  ' TOO LONG !'
            CALL EXIT_OWN(1)
          END IF
          FILE(1:)=PATH(1:I2-1)
          NFR=0
          READ (IUNIN,'(A72)') ZEILE
625       IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').EQ.0) THEN
            READ (ZEILE,6664) (DATD(IATM),IATM=1,NATMI)
            IREAD=0
            GOTO 620
          ELSE
            NFR=NFR+1
            READ (ZEILE,'(A8)') FILNAM
            WRITE (iunout,*) ' NFR =',NFR,' FILNAM = ',FILNAM
            FILE(I2:)=FILNAM
            WRITE (iunout,*) ' FILE = ',FILE
            ALLOCATE (CURFILE)
            CURFILE%RFILE = FILE
            CURFILE%NEXT => REFFILES
            REFFILES => CURFILE
            READ (IUNIN,'(A72)') ZEILE
            GOTO 625
          ENDIF
        ENDIF
      ENDIF
      READ (IUNIN,6664) (DATD(IATM),IATM=1,NATMI)
620   READ (IUNIN,6664) (DMLD(IMOL),IMOL=1,NMOLI)
      READ (IUNIN,6664) (DIOD(IION),IION=1,NIONI)
      READ (IUNIN,6664) (DPLD(IPLS),IPLS=1,NPLSI)
      IF (NPHOTI > 0)
     .  READ (IUNIN,6664) (DPHD(IPHOT),IPHOT=1,NPHOTI)

      IF (ASSOCIATED(REFFILES)) THEN
        NHD6 = NFR
      ELSE
        NHD6 = 8
      END IF
      NH0=NHD1*NHD2*NHD6
      NH1=NH0*NHD3
      NH2=NH1*NHD4
      NH3=NH2*NHD5
      CALL ALLOC_CREF
      CALL ALLOC_CREFMOD
      NFLR = 0
      DO WHILE (ASSOCIATED(REFFILES))
        NFLR = NFLR + 1
        REFFIL(NFLR) = REFFILES%RFILE
        CURFILE => REFFILES
        REFFILES => REFFILES%NEXT
        DEALLOCATE(CURFILE)
      END DO

      READ (IUNIN,6664) ERMIN,ERCUT,RPROB0,RINTEG(1),EINTEG(1),AINTEG(1)
      DO
        IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3) .EQ. '***') THEN
          IREAD = 1
          EXIT
        END IF
        IF (ZEILE(1:8) /= 'SURFMOD_') CYCLE
        ALLOCATE (REFCUR)
        ALLOCATE (REFCUR%JSRS(nspz))
        ALLOCATE (REFCUR%JSRC(nspz))
        ALLOCATE (REFCUR%TRANSPR(nspz,2))
        ALLOCATE (REFCUR%RCYCFR(nspz))
        ALLOCATE (REFCUR%RCYCTR(nspz))
        ALLOCATE (REFCUR%RCPRMR(nspz))
        ALLOCATE (REFCUR%EXPPLR(nspz))
        ALLOCATE (REFCUR%EXPELR(nspz))
        ALLOCATE (REFCUR%EXPILR(nspz))
        ALLOCATE (REFCUR%RCYCSR(nspz))
        ALLOCATE (REFCUR%RCYCCR(nspz))
        ALLOCATE (REFCUR%STPRMR(nspz))
        REFCUR%REFNAME = TRIM(ADJUSTL(ZEILE(9:)))
        IREAD=0
        READ (IUNIN,6666) REFCUR%JLREF,REFCUR%JLSPT,
     .                    REFCUR%JSRS(1),REFCUR%JSRC(1)
        READ (IUNIN,6664) REFCUR%ZNMLR,REFCUR%EWALLR,REFCUR%EWBINR,
     .                    REFCUR%TRANSPR(1,1:2),REFCUR%FSHEATR
        READ (IUNIN,6664) REFCUR%RCYCFR(1),REFCUR%RCYCTR(1),
     .                    REFCUR%RCPRMR(1),REFCUR%EXPPLR(1),
     .                    REFCUR%EXPELR(1),REFCUR%EXPILR(1)
        DO I=2,NSPZ
          REFCUR%JSRS(I) = REFCUR%JSRS(1)
          REFCUR%JSRC(I) = REFCUR%JSRC(1)
          REFCUR%TRANSPR(I,1:2)=REFCUR%TRANSPR(1,1:2)
          REFCUR%RCYCFR(I) = REFCUR%RCYCFR(1)
          REFCUR%RCYCTR(I) = REFCUR%RCYCTR(1)
          REFCUR%RCPRMR(I) = REFCUR%RCPRMR(1)
          REFCUR%EXPPLR(I) = REFCUR%EXPPLR(1)
          REFCUR%EXPELR(I) = REFCUR%EXPELR(1)
          REFCUR%EXPILR(I) = REFCUR%EXPILR(1)
        END DO
C
C  DEFAULT
        REFCUR%RCYCSR=RECYCS(1,0)
        REFCUR%RCYCCR=RECYCC(1,0)
        REFCUR%STPRMR=SPTPRM(1,0)
C
        READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        ideflt_sput=-1
        ideflt_spez=-1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
        do while (ZEILE(1:1).NE.'*'.AND.ZEILE(1:8).NE.'SURFMOD_')
          varname = zeile(1:8)
          call uppercase (varname)
          spcname = zeile(10:17)
          call uppercase (spcname)
          ispz=-1
          ico=0
          do is=1,nspz
            if (texts(is) == spcname) then
              if (ico == 0) then
                ispz = is
                ico = 1
              else
                call leer(2)
                write (iunout,*) ' WARNING !! '
                write (iunout,*) 
     .             ' ambigous species names found in surface ',
     .             'model "',REFCUR%REFNAME,'"'
                write (iunout,*) varname,'is modified for species ',ispz
              end if
            end if
          end do     ! end of "is"-loop
          ideflt_spez=1   !tentatively assume: species card
          select case (varname)
          case ('ISRS')
            if (ispz > 0)
     .          read (zeile(18:),'(I6)') REFCUR%JSRS(ispz)
          case ('ISRC')
            if (ispz > 0)
     .          read (zeile(18:),'(I6)') REFCUR%JSRC(ispz)
          case ('TRANSP1')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%TRANSPR(ispz,1)
          case ('TRANSP2')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%TRANSPR(ispz,2)
          case ('RECYCF')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCFR(ispz)
          case ('RECYCT')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCTR(ispz)
          case ('RECPRM')
            if (ispz >0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCPRMR(ispz)
          case ('EXPPL')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%EXPPLR(ispz)
          case ('EXPEL')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%EXPELR(ispz)
          case ('EXPIL')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%EXPILR(ispz)
          case ('RECYCS')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCSR(ispz)
          case ('RECYCC')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCCR(ispz)
          case ('SPTPRM')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%STPRMR(ispz)

          case default
c  not a species card, hence: a sputer model card
            if (ispz < 0) then
              READ (ZEILE,6664) REFCUR%RCYCSR(1),REFCUR%RCYCCR(1),
     .                          REFCUR%STPRMR(1)
              DO I=2,NSPZ
                REFCUR%RCYCSR(I) = REFCUR%RCYCSR(1)
                REFCUR%RCYCCR(I) = REFCUR%RCYCCR(1)
                REFCUR%STPRMR(I) = REFCUR%STPRMR(1)
              ENDDO
              ideflt_sput=1
              ideflt_spez=-1
            end if
          end select
          if ((ideflt_spez > 0) .and. (ispz < 0)) then
            write (iunout,*) 
     .        ' wrong card in species dep. reflection model'
            write (iunout,*) zeile
          end if
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
C
        END do    ! end of do-while loop (search for * or for SURFMOD_)
C all lines of this particular SURFMOD_... are read now
C
        IF (ideflt_sput.le.0.and.REFCUR%JLSPT.NE.0) THEN
          WRITE (iunout,*) 'WARNING: SPUTTERING FOR MODEL ',
     .                      REFCUR%REFNAME
          WRITE (iunout,*) 'BUT NO PARAMETERS RECYCS, RECYCC ARE READ '
          WRITE (iunout,*) 'DEFAULT MODEL: "NO SPUTTERING" IS USED. '
          WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
          REFCUR%JLSPT=0
        ENDIF
        L=LEN_TRIM(REFCUR%REFNAME)
        WRITE (iunout,*) 'SURFACE MODEL "',REFCUR%REFNAME(1:L),
     .                   '" DEFINED'
        REFCUR%NEXT => REFLIST
        REFLIST => REFCUR
      END DO
C
C  READ DATA FOR PRIMARY SOURCE  700--799
C
700   CONTINUE
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL MASAGE ('*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA   ')
C
710   READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1) .EQ. '*') GOTO 710
      READ (ZEILE,6666) NSTRAI
      IREAD=0
      WRITE (iunout,*) '       NSTRAI= ',NSTRAI
      CALL LEER(1)
!pb      READ (IUNIN,6666) (INDSRC(ISTRA),ISTRA=1,NSTRAI)
      READ (IUNIN,6666) (INDSRC(IST),IST=1,NSTRAI)
      READ (IUNIN,6664) ALLOC
      DO 712 ISTRA=1,NSTRAI
        IF (INDSRC(ISTRA).EQ.6) GOTO 712
        I=ISTRA
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') TXTSOU(ISTRA)
        ELSEIF (IREAD.EQ.1) THEN
          READ (ZEILE,'(A72)') TXTSOU(ISTRA)
          IREAD=0
        ENDIF
713     READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') GOTO 713
        READ (ZEILE,6665) NLAVRP(ISTRA),NLAVRT(I),NLSYMP(I),NLSYMT(I),
     .                    NLRAY(ISTRA)
        IREAD=0
        READ (IUNIN,6666) NPTS(ISTRA),NINITL(I),NEMODS(I),NAMODS(I),
     .                    NRAYEN(ISTRA)
csw if npts < 0 --> npts = infty
        if(npts(istra).lt.0) npts(istra) = huge(1)-1
        READ (IUNIN,66662) FLUX(ISTRA),SCALV(ISTRA),IVLSF(ISTRA),
     .                    ISCLS(ISTRA),ISCLT(ISTRA),ISCL1(ISTRA),
     .                    ISCL2(ISTRA),ISCL3(ISTRA),ISCLB(ISTRA),
     .                    ISCLA(ISTRA)
        READ (IUNIN,6665) NLATM(ISTRA),NLMOL(ISTRA),NLION(I),NLPLS(I)
        READ (IUNIN,6666) NSPEZ(ISTRA)
        READ (IUNIN,6665) NLPNT(ISTRA),NLLNE(ISTRA),
     .                    NLSRF(ISTRA),NLVOL(ISTRA),NLCNS(ISTRA)
C  SAME FOR POINT-, LINE-, SURFACE- AND VOLUME SOURCES
        READ (IUNIN,6666) NSRFSI(ISTRA)
        IF (NSRFSI(ISTRA).GT.NSRFS)
     .    CALL MASPRM('NSRFS',5,NSRFS,'NSRFSI(I)',9,NSRFSI(I),IERROR)
        DO 715 J=1,NSRFSI(ISTRA)
          READ (IUNIN,6666) INUM,INDIM(J,I),INSOR(J,I),
     .                      INGRDA(J,I,1),INGRDE(J,I,1),
     .                      INGRDA(J,I,2),INGRDE(J,I,2),
     .                      INGRDA(J,I,3),INGRDE(J,I,3)
          READ (IUNIN,6664) SORWGT(J,I),SORLIM(J,I),
     .                      SORIND(J,I),SOREXP(J,I),SORIFL(J,I)
          READ (IUNIN,6666) NRSOR(J,I),NPSOR(J,I),NTSOR(J,I),
     .                      NBSOR(J,I),NASOR(J,I),NISOR(J,I),ISTOR(J,I)
          READ (IUNIN,6664) SORAD1(J,I),SORAD2(J,I),SORAD3(J,I),
     .                      SORAD4(J,I),SORAD5(J,I),SORAD6(J,I)
715     CONTINUE
C  VELOCITY SPACE DISTRIBUTION
        READ (IUNIN,6664) SORENI(ISTRA),SORENE(I),SORVDX(I),SORVDY(I),
     .                    SORVDZ(ISTRA)
        READ (IUNIN,6664) SORCOS(I),SORMAX(I),SORCTX(I),SORCTY(I),
     .                    SORCTZ(ISTRA),RAYFRAC(ISTRA)
C
712   CONTINUE
C
      READ (IUNIN,*)
C     READ ADDITIONAL DATA FOR SOME SPECIFIC ZONES
C
800   CONTINUE
C
      CALL MASAGE ('*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES ')
810   READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1) .EQ. '*') GOTO 810
      READ (ZEILE,6666) NZADD
      IREAD=0
      WRITE (iunout,*) '       NZADD= ',NZADD
      CALL LEER(1)
      NULLIFY(TEMPLIST)
      NULLIFY(DENLIST)
      NULLIFY(VELLIST)
      NULLIFY(VOLLIST)
      DO 811 I=1,NZADD
        IF (IREAD.EQ.0) READ (IUNIN,*)
        IREAD=0
        READ (IUNIN,6666) INI,INE
        IF (INI.GT.NRAD.OR.INI.LE.0) GOTO 998
        IF (INE.GT.NRAD) GOTO 998
        IF (INE.LE.0) INE=INI
812     READ (IUNIN,'(A72)') ZEILE
C  IGJUM3 FLAG
        IF (ZEILE(1:3).EQ.'CH3') THEN
          INILGJ=MAX(1,MIN(NOPTIM,INI))
          INELGJ=MIN(NOPTIM,INE)
          DO 821 IN=INILGJ,INELGJ
          IF (NLIMPB >= NLIMPS) THEN
            CALL DEKEY (ZEILE(4:72),IGJUM3,0,NOPTIM,IN,NLIMPS)
          ELSE
            CALL DEKEYB (ZEILE(4:72),IGJUM3,0,NOPTIM,IN,NLIMPB,NBITS)
          END IF
821         CONTINUE
          GOTO 812
C  TEMPERATURE
        ELSEIF (ZEILE(1:1).EQ.'T') THEN
          DO 822 IN=INI,INE
            ALLOCATE(TEMPCUR)
            READ (ZEILE(2:72),66664) TEMPCUR%IDION,TEMPCUR%TE,
     .                               TEMPCUR%TI
            TEMPCUR%IN = IN
            TEMPCUR%NEXT => TEMPLIST
            TEMPLIST => TEMPCUR
822       CONTINUE
          GOTO 812
C  DENSITY
        ELSEIF (ZEILE(1:1).EQ.'D') THEN
          DO 823 IN=INI,INE
            ALLOCATE(DENCUR)
            READ (ZEILE(2:72),66664) DENCUR%IDION,DENCUR%DI
            DENCUR%IN = IN
            DENCUR%NEXT => DENLIST
            DENLIST => DENCUR
823       CONTINUE
          GOTO 812
C  VELOCITY (CM/SEC OR MACH)
        ELSEIF ((ZEILE(1:1).EQ.'V'.AND.ZEILE(2:2).NE.'L')
     .       .OR.ZEILE(1:1).EQ.'M') THEN
          DO 824 IN=INI,INE
            ALLOCATE(VELCUR)
            READ (ZEILE(2:72),66664) VELCUR%IDION,VELCUR%VX,VELCUR%VY,
     .                               VELCUR%VZ
            IF (ZEILE(1:1).EQ.'M') THEN
              VELCUR%IZ=1
            ELSEIF (ZEILE(1:1).EQ.'V') THEN
              VELCUR%IZ=-1
            ENDIF
            VELCUR%IN = IN
            VELCUR%NEXT => VELLIST
            VELLIST => VELCUR
824       CONTINUE
          GOTO 812
C  VOLUME
        ELSEIF (ZEILE(1:2).EQ.'VL') THEN
          DO 825 IN=INI,INE
            ALLOCATE(VOLCUR)
            READ (ZEILE(3:72),6664) VOLCUR%VOL
            VOLCUR%IN = IN
            VOLCUR%NEXT => VOLLIST
            VOLLIST => VOLCUR
825       CONTINUE
          GOTO 812
        ELSEIF (ZEILE(1:1).EQ.'*') THEN
          IREAD=1
          GOTO 811
        ELSE
          GOTO 998
        ENDIF
811   CONTINUE

C
C  READ DATA FOR STATISTICS AND NONANALOG MODEL, 900--999
C
900   CONTINUE
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      IREAD=0
      CALL MASAGE ('*** 9. DATA FOR STATISTIC AND NONANALOG MODEL   ')
C
910   READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 910
C  DATA FOR CONDITIONAL EXPECTATION ESTIMATOR
      READ (ZEILE,6665) (NLPRCA(J),J=1,NATM),(NLPRCM(J),J=1,NMOL),
     .                  (NLPRCI(J),J=1,NION),(NLPRCPH(J),J=1,NPHOT)
      READ (IUNIN,6666) NPRCSF
      NPRCSF=MIN0(NLIMPS,NPRCSF)
      IPRCSF=1
911   CONTINUE
      IF (IPRCSF.LE.NPRCSF) THEN
        READ (IUNIN,6666) (IPRSF(J),J=1,12)
        DO J=1,12
        IF (IPRSF(J).GT.0.AND.IPRSF(J).LE.NLIMPS)
     .     NLPRCS(IPRSF(J))=.TRUE.
        ENDDO
        IPRCSF=IPRCSF+12
        GOTO 911
      ENDIF
C  DATA FOR SPLITTING AND RUSSIAN ROULETTE
      READ (IUNIN,6666) MAXLEV,MAXRAD,MAXPOL,MAXTOR,MAXADD
      MXL=15
      IF (MAXLEV.GT.MXL)
     .    CALL MASPRM('MXL',3,MXL,'MAXLEV',6,MAXLEV,IERROR)
      IF (ABS(MAXRAD).GT.N1ST)
     .    CALL MASPRM('N1ST',4,N1ST,'MAXRAD',6,MAXRAD,IERROR)
      IF (MAXPOL.GT.N2ND)
     .    CALL MASPRM('N2ND',4,N2ND,'MAXPOL',6,MAXPOL,IERROR)
      IF (MAXTOR.GT.N3RD)
     .    CALL MASPRM('N3RD',4,N3RD,'MAXTOR',6,MAXTOR,IERROR)
      IF (MAXADD.GT.NLIM)
     .    CALL MASPRM('NLIM',4,NLIM,'MAXADD',6,MAXADD,IERROR)
      DO IN=1,MAXRAD
        READ (IUNIN,66665) ID,NSSPL(IN),PRMSPL(IN)
      ENDDO
      DO IN=1,MAXPOL
        READ (IUNIN,66665) ID,NSSPL(N1ST+IN),PRMSPL(N1ST+IN)
      ENDDO
      DO IN=1,MAXTOR
        READ (IUNIN,66665) ID,NSSPL(N1ST+N2ND+IN),PRMSPL(N1ST+N2ND+IN)
      ENDDO
      DO IN=1,MAXADD
        READ (IUNIN,66665) ID,NSSPL(N1ST+N2ND+N3RD+IN),
     .                        PRMSPL(N1ST+N2ND+N3RD+IN)
      ENDDO
C  DATA FOR BIAS SAMPLING
      READ (IUNIN,6664) WMINV,WMINS,WMINC,WMINL
      WMINV=MAX(WMINV,EPS60)
      WMINS=MAX(WMINS,EPS60)
      WMINC=MAX(WMINC,EPS60)
      READ (IUNIN,6664) SPLPAR
C  DATA FOR STANDARD DEVIATION
      READ (IUNIN,*)
      CALL LEER(1)
      WRITE (iunout,*) '       CARDS FOR STANDARD DEVIATION '
      READ (IUNIN,6666) NSIGVI,NSIGSI,NSIGCI,NSIGI_BGK,NSIGI_COP,
     .                  NSIGI_SPC
      WRITE (iunout,*) '       NSIGVI,NSIGSI,NSIGCI= ',
     .                         NSIGVI,NSIGSI,NSIGCI
      WRITE (iunout,*) '       NSIGI_BGK,NSIGI_COP = ',
     .                         NSIGI_BGK,NSIGI_COP
      WRITE (iunout,*) '       NSIGI_SPC           = ',NSIGI_SPC
      CALL LEER(1)
      DO 913 J=1,NSIGVI
        READ (IUNIN,6666) IGH(J),IIH(J)
913   CONTINUE
      DO 914 J=1,NSIGSI
        READ (IUNIN,6666) IGHW(J),IIHW(J)
914   CONTINUE
      DO 915 J=1,NSIGCI
        READ (IUNIN,6666) IGHC(1,J),IIHC(1,J),IGHC(2,J),IIHC(2,J)
915   CONTINUE
C
C   READ DATA FOR ADDITIONAL AND SURFACE AVERAGED TALLIES
      READ (IUNIN,*)
      CALL MASAGE ('*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION   ')
      CALL MASAGE ('        ESTIMATORS  AND ALGEBRAIC EXPRESSIONS    ')
C
1010  READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 1010
      READ (ZEILE,6666) NADVI,NCLVI,NALVI,NADSI,NALSI,NADSPC
      WRITE (iunout,*) '        NADVI,NCLVI,NALVI= ',
     .                     NADVI,NCLVI,NALVI
      WRITE (iunout,*) '        NADSI,NALSI,NADSPC= ',
     .                     NADSI,NALSI,NADSPC
      CALL LEER(1)
C
C
1000  CONTINUE
C
      READ (IUNIN,*)
      CALL MASAGE('*** 10A. DATA FOR ADDITIONAL TALLIES           ')
      ALLOCATE (TXTTLA(NADVI))
      ALLOCATE (TXTSCA(NADVI))
      ALLOCATE (TXTUTA(NADVI))
      DO 1020 J=1,NADVI
1021    READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1021
        READ (ZEILE,6666) IADVE(J),IADVS(J),IADVT(J),IADRC(J)
        READ (IUNIN,'(A72)') TXTTLA(J)
        READ (IUNIN,'(2A24)') TXTSCA(J),TXTUTA(J)
1020  CONTINUE
      READ (IUNIN,*)
      CALL MASAGE('*** 10B. DATA FOR COLLISION ESTIMATORS         ')
      ALLOCATE (TXTTLC(NCLVI))
      ALLOCATE (TXTSCC(NCLVI))
      ALLOCATE (TXTUTC(NCLVI))
      DO 1030 J=1,NCLVI
1031    READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1031
        READ (ZEILE,6666) ICLVE(J),ICLVS(J),ICLVT(J),ICLRC(J)
        READ (IUNIN,'(A72)') TXTTLC(J)
        READ (IUNIN,'(2A24)') TXTSCC(J),TXTUTC(J)
1030  CONTINUE
      READ (IUNIN,*)
      CALL MASAGE('*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS     ')
      ALLOCATE (TXTTLR(NALVI))
      ALLOCATE (TXTSCR(NALVI))
      ALLOCATE (TXTUTR(NALVI))
      DO 1040 J=1,NALVI
1041    READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1041
        READ (ZEILE,'(A72)') CHRTAL(J)
        READ (IUNIN,'(A72)') TXTTLR(J)
        READ (IUNIN,'(2A24)') TXTSCR(J),TXTUTR(J)
1040  CONTINUE
C
      READ (IUNIN,*)
      CALL MASAGE('*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES   ')
      DO 1050 J=1,NADSI
1051    READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1051
        READ (ZEILE,6666) IADSE(J),IADSS(J),IADST(J),IADSC(J)
        READ (IUNIN,'(A72)') TXTTLW(J,NTLSA)
        READ (IUNIN,'(2A24)') TXTSPW(J,NTLSA),TXTUNW(J,NTLSA)
1050  CONTINUE
C
      READ (IUNIN,*)
      CALL MASAGE('*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES ')
      DO 1060 J=1,NALSI
1061    READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1061
        READ (ZEILE,'(A72)') CHRTLS(J)
        READ (IUNIN,'(A72)') TXTTLW(J,NTLSR)
        READ (IUNIN,'(2A24)') TXTSPW(J,NTLSR),TXTUNW(J,NTLSR)
1060  CONTINUE
C
      READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:3) == '** ') THEN
        CALL MASAGE('*** 10F. DATA FOR SPECTRA                   ')
        IF (NADSPC > 0) THEN
          ALLOCATE(ESTIML(NADSPC))
          IF (NSMSTRA > 0) ALLOCATE(SMESTL(NADSPC))
        END IF
        DO J=1,NADSPC
          DO
            READ (IUNIN,'(A72)') ZEILE
            IF (ZEILE(1:1) .NE. '*') EXIT
          END DO
          READ (ZEILE,'(12I6)') ISPSRF, IPTYP, IPSP, ISPTYP, NSPS,
     .                          ISRFCLL, IDIREC
          READ (IUNIN,'(6E12.4)') SPCMN, SPCMX, SPC_SHIFT,
     .                             SPCPLT_X,SPCPLT_Y,SPCPLT_SAME
          SPCVX = 0._DP
          SPCVY = 0._DP
          SPCVZ = 0._DP
          IF (IDIREC /= 0) THEN
            READ (IUNIN,'(6E12.4)') SPCVX, SPCVY, SPCVZ
            VNORM = SQRT(SPCVX**2+SPCVY**2+SPCVZ**2)+EPS60
            SPCVX = SPCVX / VNORM
            SPCVY = SPCVY / VNORM
            SPCVZ = SPCVZ / VNORM
            IF (ISRFCLL == 0) THEN
              WRITE (IUNOUT,*) ' SPECTRUM NUMBER ',J
              WRITE (IUNOUT,*) ' DEFINITION OF LINE OF SIGHT FOR',
     .                         ' SURFACE SPECTRUM IS NOT FORSEEN '
              WRITE (IUNOUT,*) ' SPECIFICATION OF DIRECTION IS',
     .                         ' IGNORED '
              IDIREC = 0
            END IF
          END IF 
          IF (ISPSRF > 0) THEN
            IF ((ISRFCLL == 0) .AND. (ISPSRF > NLIMI)) THEN
              WRITE (iunout,*) 
     .          ' SURFACE INDEX FOR SPECTRUM OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' SURFACE NUMBER = ',ISPSRF
              IERROR = IERROR + 1
            END IF
            IF (((ISRFCLL == 1) .AND. (ISPSRF > NRTAL)) .OR. 
     .          ((ISRFCLL == 2) .AND. (ISPSRF > NRAD))) THEN
              WRITE (iunout,*) 
     .          ' CELL INDEX FOR SPECTRUM OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' CELL NUMBER = ',ISPSRF
              IERROR = IERROR + 1
            END IF
          ELSEIF (ISPSRF < 0) THEN
            ISPSRF = ABS(ISPSRF)
            IF (ISPSRF > NSTSI) THEN
              WRITE (iunout,*) 
     .          ' SURFACE INDEX FOR SPECTRUM OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' SURFACE NUMBER = ',ISPSRF
              IERROR = IERROR + 1
            END IF
            ISPSRF = NLIM+ISPSRF
          ELSE
            WRITE (iunout,*) ' SURFACE INDEX 0 NOT FORSEEN '
            WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
            IERROR = IERROR + 1
          END IF

          IF ((IPTYP < 0) .OR. (IPTYP > 4)) THEN
            WRITE (iunout,*) ' PARTICLE TYPE ',IPTYP,' NOT FORSEEN '
            WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
            IERROR = IERROR + 1
          ELSE
            IF (((IPTYP == 0).AND.((IPSP < 0).OR.(IPSP > NPHOTI))) .OR.
     .          ((IPTYP == 1).AND.((IPSP < 0).OR.(IPSP > NATMI))) .OR.
     .          ((IPTYP == 2).AND.((IPSP < 0).OR.(IPSP > NMOLI))) .OR.
     .          ((IPTYP == 3).AND.((IPSP < 0).OR.(IPSP > NIONI))) .OR.
     .          ((IPTYP == 4).AND.((IPSP < 0).OR.(IPSP > NPLSI)))) THEN
              WRITE (iunout,*) ' PARTICLE SPECIES INDEX OUT OF BOUNDS '
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' SPECIES NUMBER = ',IPSP
              IERROR = IERROR + 1
            END IF
          END IF

          IF ((ISRFCLL == 0) .AND. (ISPTYP > 2)) THEN
            WRITE (IUNOUT,*) ' WRONG TYPE OF SPECTRUM SPECIFIED '
            WRITE (IUNOUT,*) ' SPECTRUM NUMBER = ',J
            WRITE (IUNOUT,*) ' SPECTRUM TYPE = ',ISPTYP
          END IF

          ALLOCATE(ESPEC)
          ESPEC%ISPCSRF = ISPSRF
          ESPEC%IPRTYP = IPTYP
          ESPEC%IPRSP = IPSP
          ESPEC%ISPCTYP = ISPTYP
          ESPEC%NSPC = NSPS
          ESPEC%IMETSP = 0
          ESPEC%ISRFCLL = ISRFCLL
          ESPEC%IDIREC = IDIREC
          ESPEC%SPCMIN=SPCMN
          ESPEC%SPCMAX=SPCMX
          ESPEC%ESP_00=SPC_SHIFT
          ESPEC%SPC_XPLT=SPCPLT_X
          ESPEC%SPC_YPLT=SPCPLT_Y
          ESPEC%SPC_SAME=SPCPLT_SAME
          ESPEC%SPCVX=SPCVX
          ESPEC%SPCVY=SPCVY
          ESPEC%SPCVZ=SPCVZ
          ESPEC%ESP_MIN=1.E30
          ESPEC%ESP_MAX=-1.E30
          ESPEC%SPCDEL=(SPCMX-SPCMN)/REAL(NSPS,DP)
          ESPEC%SPCDELI=1._DP/ESPEC%SPCDEL
          ALLOCATE(ESPEC%SPC(0:NSPS+1))
!          IF (NSIGI_SPC > 0) THEN
            ALLOCATE(ESPEC%SDV(0:NSPS+1))
            ALLOCATE(ESPEC%SGM(0:NSPS+1))
!          END IF
          ESPEC%SPC(0:NSPS+1) = 0._DP
          IF (NSMSTRA > 0) THEN
            ALLOCATE(SSPEC)
            ALLOCATE(SSPEC%SPC(0:NSPS+1))
!            IF (NSIGI_SPC > 0) THEN
              ALLOCATE(SSPEC%SDV(0:NSPS+1))
              ALLOCATE(SSPEC%SGM(0:NSPS+1))
!            END IF
            SSPEC = ESPEC
            SMESTL(J)%PSPC => SSPEC
          END IF
          ESTIML(J)%PSPC => ESPEC
        END DO
        IREAD=0
      END IF
C
C   READ DATA FOR NUMERICAL AND GRAPHICAL OUTPUT 1100--1199
C
1100  CONTINUE
C
      IF (IREAD == 0) READ (IUNIN,*)
      CALL LEER(1)
      CALL MASAGE ('*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT  ')
c
c  search for input block 11a
1110  READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 1110
      READ (ZEILE,6665) TRCPLT,TRCHST,TRCNAL,TRCMOD,TRCSIG,
     .                  TRCGRD,TRCSUR,TRCREF,TRCFLE,TRCAMD,
     .                  TRCINT,TRCLST,TRCSOU,TRCREC,TRCTIM,
     .                  TRCBLA,TRCBLM,TRCBLI,TRCBLP,TRCBLE,
     .                  TRCBLPH,TRCTAL
      READ (IUNIN,6665) (TRCSRC(J),J=0,NSTRA)
C
      READ (IUNIN,6666) NVOLPR, NSPCPR
      NVOLPR=MIN(NVOLPR,100)
      WRITE (iunout,*) '        NVOLPR= ',NVOLPR
      CALL LEER(1)
      DO 1120 J=1,NVOLPR
        READ (IUNIN,6666) NTLV,NFLGV,NSPZV1,NSPZV2,NTLVF
        IF (NTLV.LT.-NTALI.OR.NTLV.GT.NTALV) GOTO 990
        IF ((NTLVF.NE.0).AND.(NTLVF.LT.70)) GOTO 995
        NPRTLV(J)  =NTLV
        NFLAGV(J)  =NFLGV
        NSPEZV(J,1)=NSPZV1
        NSPEZV(J,2)=NSPZV2
        NTLVFL(J)  =NTLVF
1120  CONTINUE
C
      READ (IUNIN,6666) NSURPR
      IF(NSURPR > NLIMPS) THEN
         WRITE (iunout,*) 'NSURPR > NLIMPS, EXIT.'
         WRITE (iunout,*) ' NSURPR, NLIMPS ', NSURPR, NLIMPS
         CALL EXIT_OWN(1)
      ENDIF
      NSURPR=MIN(NSURPR,NLIMPS)
      WRITE (iunout,*) '        NSURPR= ',NSURPR
      CALL LEER(1)
      DO 1130 J=1,NSURPR
        READ (IUNIN,6666) NSRF,NTLS,NFLGS,NSPZS1,NSPZS2,NTLSF
        IF (NSRF.LT.0) NSRF=NLIM+IABS(NSRF)
        IF (NSRF.EQ.0.AND.NLIM+NSTSI.LT.NLIMPS) NSRF=NLIM+NSTSI+1
        IF (NSRF.LE.0.OR.NSRF.GT.NLIMPS) GOTO 991
        NPRSRF(J)=NSRF
        NPRTLS(J)=NTLS
        NFLAGS(J)=NFLGS
        NSPEZS(J,1)=NSPZS1
        NSPEZS(J,2)=NSPZS2
        NTLSFL(J)=NTLSF
1130  CONTINUE

c  overrule default switching on/off of volume averaged tallies

      READ (IUNIN,'(A72)') ZEILE
      CALL UPPERCASE(ZEILE)
      IREAD=1
      IF ((ZEILE(1:1) .NE.'*') .AND. (SCAN(ZEILE,'FT') == 0)) THEN

c  there are ntlvout volume tallies to be dealt with
c  and ntlsout surface tallies

        IREAD=0
        READ (ZEILE,6666) NTLVOUT
        ITLVOUT=0
        DO WHILE ((NTLVOUT > 0) .AND. (ITLVOUT < NTLVOUT))
          READ (IUNIN,6666) (NUMTAL(J),J=1,12)
          DO J=1,12
            IF ((NUMTAL(J) /= 0) .AND. (ABS(NUMTAL(J)) <= NTALV)) THEN
              IF (NUMTAL(J) < 0) THEN
                LMISTALV(ABS(NUMTAL(J))) = .TRUE.  ! SWITCHED OFF
              ELSE
                LMISTALV(NUMTAL(J)) = .FALSE.      ! SWITCHED ON
              END IF
            END IF
          END DO
          ITLVOUT=ITLVOUT+12
        END DO

c  overrule default switching on/off of surface averaged tallies

        READ (IUNIN,6666) NTLSOUT
        ITLSOUT=0
        DO WHILE ((NTLSOUT > 0) .AND. (ITLSOUT < NTLSOUT))
          READ (IUNIN,6666) (NUMTAL(J),J=1,12)
          DO J=1,12
            IF ((NUMTAL(J) /= 0) .AND. (ABS(NUMTAL(J)) <= NTALS)) THEN
              IF (NUMTAL(J) < 0) THEN
                LMISTALS(ABS(NUMTAL(J))) = .TRUE.  ! SWITCHED OFF
              ELSE
                LMISTALS(NUMTAL(J)) = .FALSE.      ! SWITCHED ON
              END IF
            END IF
          END DO
          ITLSOUT=ITLSOUT+12
        END DO
      END IF


C  search for input block 11b

1131  IF (IREAD == 0) READ (IUNIN,'(A72)') ZEILE
      IREAD=0
      IF (ZEILE(1:1) .EQ. '*') GOTO 1131
C  2D GEOMETRY PLOT
      READ (ZEILE,6665) PL1ST,PL2ND,PL3RD,PLADD,PLHST,
     .                  PLCUT(1),PLCUT(2),PLCUT(3),PLBOX,PLSTOR,
     .                  PLNUMV,PLNUMS,PLARR,LRPSCUT
      READ (IUNIN,6666) NPLINR,NPLOTR,NPLDLR,NPLINP,NPLOTP,NPLDLP,
     .                  NPLINT,NPLOTT,NPLDLT
C  3D GEOMETRY PLOT
      DO 1140 J=1,5
        READ (IUNIN,6662) PL3A(J),TEXTLA(J),IPLTA(J),
     .                (IPLAA(J,I),IPLEA(J,I),I=1,IPLTA(J))
1140   CONTINUE
      DO 1141 J=1,3
        READ (IUNIN,6662) PL3S(J),TEXTLS(J),IPLTS(J),
     .                (IPLAS(J,I),IPLES(J,I),I=1,IPLTS(J))
1141  CONTINUE
C
      READ (IUNIN,6664) CH2MX,CH2MY,      CH2X0,CH2Y0,CH2Z0
      READ (IUNIN,6664) CH3MX,CH3MY,CH3MZ,CH3X0,CH3Y0,CH3Z0
      READ (IUNIN,6664) ANGLE1,ANGLE2,ANGLE3
C
C  PARTICLE HISTORY PLOTS IN 2D OR 3D GEOMETRY PLOTS
      READ (IUNIN,6666) I1TRC,I2TRC,(ISYPLT(J),J=1,8),ILINIE
C

c  search for input block 11c

1151  READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 1151
C  DATA FOR PLOTS OF VOLUME AVERAGED TALLIES
      READ (ZEILE,6666) NVOLPL
      WRITE (iunout,*) '        NVOLPL= ',NVOLPL
      CALL LEER(1)
C
      IF (NVOLPL.GT.NPTAL)
     .  CALL MASPRM('NPTAL',5,NPTAL,'NVOLPT',6,NVOLPL,IERROR)
C
      IF (NVOLPL.LE.0) GOTO 1175
      READ (IUNIN,6665) (PLTSRC(J),J=0,NSTRA)
      IF (LRPSCUT) THEN
        READ (IUNIN,6664) CUTPLANE(1:4)
!  TEST IF CUTPLANE IS AXIS-PARALLEL
        IF ((SUM(ABS(CUTPLANE(2:4)))-1._DP > EPS10) .OR.
     .      (CUTPLANE(2)*CUTPLANE(3) > 0._DP) .OR.
     .      (CUTPLANE(2)*CUTPLANE(4) > 0._DP) .OR.
     .      (CUTPLANE(3)*CUTPLANE(4) > 0._DP)) THEN
          WRITE (iunout,*) ' PLANE IS NOT PARALLEL TO ONE AXIS '
          WRITE (iunout,*) ' CUTTING OF TETRAHEDRONS ABANDONNED '
          WRITE (iunout,*) ' LRPSCUT SET TO .FALSE. '
          LRPSCUT = .FALSE.
        END IF
      END IF
      IPLANE=1
      LRAPS3D=.FALSE.
      LR3DCON=.FALSE.
      RAPSDEL=-HUGE(1._DP)
      DO 1150 J=1,NVOLPL
1152    READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1152
        READ (ZEILE,6666) NSP
        IF (NSP.GT.NPLT)
     .    CALL MASPRM('NPLT',4,NPLT,'NSP',3,NSP,IERROR)
        NSPTAL(J)=NSP
        READ (IUNIN,6665) PLTL2D(J),PLTL3D(J),PLTLLG(J),PLTLER(J)
        READ (IUNIN,6664) TALZMI(J),TALZMA(J),TALXMI(J),TALXMA(J),
     .                    TALYMI(J),TALYMA(J)
        IF (PLTL2D(J)) THEN
          READ (IUNIN,6665) LHIST2(J),LSMOT2(J)
          DO 1160 I=1,NSPTAL(J)
            READ (IUNIN,6666) ISPTAL(J,I),NTL,
     .                        NPLIN2(J,I),NPLOT2(J,I),NPLDL2(J,I)
            IF (NTL.LT.-NTALI.OR.NTL.GT.NTALV.OR.NTL.EQ.0) GOTO 990
            NPTALI(J,I)=NTL
1160      CONTINUE
        ENDIF
        IF (PLTL3D(J)) THEN
          READ (IUNIN,6665) LHIST3(J),LCNTR3(J),LSMOT3(J),
     .                      LRAPS3(J),LVECT3(J),LRPVC3(J),
     .                      LRPS3D,LRPSCN
          READ (IUNIN,6665) LPRAD3(J),LPPOL3(J),LPTOR3(J)
          DO 1161 I=1,NSPTAL(J)
            READ (IUNIN,6666) ISPTAL(J,I),NTL,IPROJ3(J,I),
     .                        NPLI13(J,I),NPLO13(J,I),
     .                        NPLI23(J,I),NPLO23(J,I),IPLN
            IF (NTL.LT.-NTALI.OR.NTL.GT.NTALV.OR.NTL.EQ.0) GOTO 990
            NPTALI(J,I)=NTL
1161      CONTINUE
          READ (IUNIN,6664) TALW1(J),TALW2(J),FCABS1(J),FCABS2(J),
     .                      RPSDL
          IF (FCABS1(J).LE.0.0D0) FCABS1(J)=1.0D0
          IF (FCABS2(J).LE.0.0D0) FCABS2(J)=1.0D0
          LRAPS3D=LRAPS3D.OR.LRPS3D
          LR3DCON=LR3DCON.OR.LRPSCN
          IPLANE=MAX(IPLANE,IPLN)
          RAPSDEL=MAX(RAPSDEL,RPSDL)
        ENDIF
1150  CONTINUE
      IF (NLTRA) RAPSDEL=RAPSDEL*DEGRAD
C
C  SKIP INPUT LINES, UNTIL INPUT BLOCK 12 STARTS
1175  READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:3) .EQ. '***') GOTO 1205
      GOTO 1175
C
C  READ DATA FOR DIAGNOSTIC MODULE  1200--1299
C
1200  CONTINUE
C
1205  CALL MASAGE ('*** 12. DATA FOR DIAGNOSTIC MODULE              ')
1210  READ (IUNIN,'(A72)') ZEILE
      IREAD=0
      IF (ZEILE(1:1) .EQ. '*') GOTO 1210
      READ (ZEILE,6666) NCHORI,NCHENI
      NCHOR = NCHORI
      NCHEN = NCHENI
      WRITE (iunout,*) '        NCHORI,NCHENI= ',NCHORI,NCHENI
      CALL LEER(1)
      IF (IABS(NCHENI).GT.NCHEN)
     .    CALL MASPRM('NCHEN',5,NCHEN,'IABS(NCHENI)',12,IABS(NCHENI),
     .                 IERROR)
      IF (NCHORI.LE.0) GOTO 1230
      CALL ALLOC_COMSIG
      CALL INIT_COMSIG
      DO 1220 ICHORI=1,NCHORI
        READ (IUNIN,'(A72)') TXTSIG(ICHORI)
        READ (IUNIN,6666) NCHTAL(ICHORI),NSPSCL(ICHORI),NSPNEW(ICHORI),
     .                    ISTCHR
        READ (IUNIN,6666) NSPSTR(ICHORI),NSPSPZ(ICHORI),
     .                    NSPINI(ICHORI),NSPEND(ICHORI),
     .                    NSPBLC(ICHORI),NSPADD(ICHORI)
        READ (IUNIN,6664) EMIN1(ICHORI),EMAX1(ICHORI),ESHIFT(ICHORI)
        READ (IUNIN,66664) IPIVOT(ICHORI),
     .                     XPIVOT(ICHORI),YPIVOT(ICHORI),ZPIVOT(ICHORI)
        READ (IUNIN,66664) ICHORD(ICHORI),
     .                     XCHORD(ICHORI),YCHORD(ICHORI),ZCHORD(ICHORI)
        NLSTCHR(ICHORI) = ISTCHR > 0
1220  CONTINUE
      READ (IUNIN,6665) PLCHOR,PLSPEC
1230  CONTINUE
C  SKIP READING REST OF THIS BLOCK
      READ (IUNIN,'(A72)') ZEILE
      IREAD=0
      IF (ZEILE(1:3).NE.'***') GOTO 1230
      IREAD=1
C
C  READ DATA FOR NONLINEAR MODE  1300--1399
C
1300  CONTINUE
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      IREAD=0
      CALL MASAGE ('*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION ')
C
      READ (IUNIN,6666) NPRNLI
      WRITE (iunout,*) '        NPRNLI= ',NPRNLI
      CALL LEER(1)
      IF (NLERG.AND.NPRNLI.LE.0) THEN
C  NO TIME HORIZON DEFINED, DESPITE NLERG=.TRUE.
C  THEREFORE: SET A DEFAULT TIME HORIZON HERE
        NPRNLI=100
      ENDIF
      IF (NPRNLI.LE.0) GOTO 1350
      READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1).EQ.'*') THEN
C  DATA FOR DEFAULT TIME HORIZON
        NTIME0=0
        NTIME=1
        DTIMV=1.D0
        TIME0=0.D0
        NPTST=0
        NTMSTP=1
        NSNVI=0
        GOTO 1350
      ENDIF
      READ (ZEILE,6666) NPTST,NTMSTP
      IREAD=0
      READ (IUNIN,6664) DTIMV,TIME0
C   READ DATA FOR SNAPSHOT TALLIES
      READ (IUNIN,*)
      CALL MASAGE ('*** 13A. DATA FOR SNAPSHOT TALLIES           ')
      READ (IUNIN,6666) NSNVI
      WRITE (iunout,*) '        NSNVI= ',NSNVI
      CALL LEER(1)
      ALLOCATE (TXTTLT(NSNVI))
      ALLOCATE (TXTSCT(NSNVI))
      ALLOCATE (TXTUTT(NSNVI))
      DO 1320 J=1,NSNVI
1321    READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') GOTO 1321
        READ (ZEILE,6666) ISNVE(J),ISNVS(J),ISNVT(J),ISNRC(J)
        IREAD=0
!        READ (IUNIN,'(A72)') TXTTAL(J,NTALT)
!        READ (IUNIN,'(2A24)') TXTSPC(J,NTALT),TXTUNT(J,NTALT)
        READ (IUNIN,'(A72)') TXTTLT(J)
        READ (IUNIN,'(2A24)') TXTSCT(J),TXTUTT(J)
1320  CONTINUE
C
      IF (NTIME.LE.0) THEN
        WRITE (iunout,*) 'ERROR IN INPUT: TIME DEP. MODE BUT NTIME.LE.0'
        CALL EXIT_OWN(1)
      ENDIF
1350  CONTINUE
C  SKIP READING REST OF THIS BLOCK
      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3).NE.'***') GOTO 1350
      IREAD=1
C
      IF (NTIME.GE.1) THEN
C  DEFINE ONE MORE SURFACE
        NSTSI=NSTSI+1
C  CHECK STORAGE
        CALL LEER(1)
        IF (NSTSI.GT.NSTS) THEN
          CALL MASPRM('NSTS',4,NSTS,'NSTSI',5,NSTSI,IERROR)
          CALL EXIT_OWN(1)
        ENDIF
C
C  SET DEFAULTS FOR "TIME HORIZON"
C
        TXTSFL(NLIM+NSTSI)='"TIME HORIZON"                           '
        ILIIN(NLIM+NSTSI)=2
C
C
cdr     IF (NFILEJ.EQ.2.OR.NFILEJ.EQ.3) THEN
cdr functioniert noch nicht, falls mehrere timesteps, davon nur
cdr der erste: initialisierung, die anderen: fortsetzung.
cdr denn dann wird bei der fortsetzung das stratum nicht gemacht.
cdr wg. goto 4000. angefangen: "mkcens" (make stratum for census array)
C
C  DEFINE ONE MORE STRATUM
          NSTRAI=NSTRAI+1
C  CHECK STORAGE
          IF (NSTRAI.GT.NSTRA) THEN
            CALL MASPRM('NSTRA',5,NSTRA,'NSTRAI',6,NSTRAI,IERROR)
            CALL EXIT_OWN(1)
          ENDIF
C
C  SET DEFAULTS FOR SOURCE DUE TO INITIAL CONDITION, VALID ONLY FOR
C  FIRST TIMESTEP. MODIFIED FOR LATER TIMESTEPS IN SUBR. TMSTEP
C
          TXTSOU(NSTRAI)='SOURCE DUE TO INITIAL CONDITION          '
C  SOURCE DISTRIBUTION SAMPLED FROM CENSUS ARRAYS RPARTC,IPARTC
          NLCNS(NSTRAI)=.TRUE.
          NLPNT(NSTRAI)=.FALSE.
          NLLNE(NSTRAI)=.FALSE.
          NLSRF(NSTRAI)=.FALSE.
          NLVOL(NSTRAI)=.FALSE.
C  DO NOT CALL IF2COP(NSTRAI)
          INDSRC(NSTRAI)=-1
C
          NLAVRP(NSTRAI)=.FALSE.
          NLAVRT(NSTRAI)=.FALSE.
          NLSYMP(NSTRAI)=.FALSE.
          NLSYMT(NSTRAI)=.FALSE.
          NPTS(NSTRAI)=0
          NINITL(NSTRAI)=2000*NINITL(NSTRAI-1)+1
          NEMODS(NSTRAI)=1
          NAMODS(NSTRAI)=1
          FLUX(NSTRAI)=0.
          NLATM(NSTRAI)=.FALSE.
          NLMOL(NSTRAI)=.FALSE.
          NLION(NSTRAI)=.FALSE.
          NLPLS(NSTRAI)=.FALSE.
          NSPEZ(NSTRAI)=0
          NSRFSI(NSTRAI)=0
C
          SORENI(NSTRAI)=0.
          SORENE(NSTRAI)=0.
          SORVDX(NSTRAI)=0.
          SORVDY(NSTRAI)=0.
          SORVDZ(NSTRAI)=0.
          SORCOS(NSTRAI)=0.
          SORMAX(NSTRAI)=0.
          SORCTX(NSTRAI)=0.
          SORCTY(NSTRAI)=0.
          SORCTZ(NSTRAI)=0.
C
C  NEW TIMESTEP
C
          IF (DTIMVN.LE.0.D0) THEN
            DTIMVN=DTIMV
C         ELSE
C           DTIMVN=DTIMVN
          ENDIF
C
C  OLD TIMESTEP
C
C  READ INITIAL POPULATION FROM PREVIOUS RUN, OVERWRITE DEFAULTS
C
        IPRNL=0
        IF (NFILEJ.EQ.2.OR.NFILEJ.EQ.3) THEN
          CALL RSNAP
          DTIMVO=DTIMV
C
          WRITE (iunout,*) 'INITIAL POPULATION FOR FIRST TIMESTEP'
          WRITE (iunout,*) 'READ FROM FILE FT 15 '
          WRITE (iunout,*) 'PARTICLES AND FLUX STORED FOR INITIAL '
          WRITE (iunout,*) 'DISTRIBUTION IN PREVIOUS RUN '
          CALL MASJ1('IPRNL   ',IPRNL)
          CALL MASR1('FLUX    ',FLUX(NSTRAI))
C
          IF (DTIMVN.NE.DTIMVO) THEN
            FLUX(NSTRAI)=FLUX(NSTRAI)*DTIMVO/DTIMVN
C
            WRITE (iunout,*) 'FLUX IS RESCALED BY DTIMV_OLD/DTIMV_NEW '
            CALL MASR1('FLUX    ',FLUX(NSTRAI))
            CALL LEER(1)
          ENDIF
C
          CALL LEER(2)
          IF (TIME0.GT.0.) THEN
            DO I=1,IPRNL
              RPARTC(I,10)=TIME0
            ENDDO
            WRITE (iunout,*) 'PARTICLE CLOCK RESET TO TIME0'
            WRITE (iunout,*) 'FIRST TIMESTEP RUNS FROM TIM1 TO TIM2:  '
            CALL MASR2('TIM1, TIM2      ',TIME0,TIME0+DTIMV)
            CALL LEER(2)
          ENDIF
C
        ENDIF
C
        DTIMV=DTIMVN
C
        IF (NPTST.EQ.0) THEN
          NPTS(NSTRAI)=IPRNL
          NLMOVIE=.FALSE.
        ELSEIF (NPTST.GT.0) THEN
          NPTS(NSTRAI)=NPTST
          NLMOVIE=.FALSE.
        ELSEIF (NPTST.LT.0) THEN
          NPTS(NSTRAI)=IPRNL
          NLMOVIE=.TRUE.
        ENDIF
C
        IF (NPTS(NSTRAI).GT.0.AND.FLUX(NSTRAI).GT.0) THEN
          NSRFSI(NSTRAI)=1
          SORWGT(1,NSTRAI)=1.D0
        ENDIF
C
      ENDIF
1399  CONTINUE
C
1500  IF (IERROR.GT.0) THEN
        WRITE (iunout,*) IERROR,' INPUT OR PARAMETER ERRORS DETECTED '
        WRITE (iunout,*) 
     .    ' SEE THE ERRORMESSAGES LISTED ABOVE AND CORRECT '
        WRITE (iunout,*) ' THE ERRORS BEFORE RE-EXECUTION '
        CALL EXIT_OWN(1)
      ENDIF
C
      CALL PAGE

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time nach einlesen ',tpb2-tpb1
!pb      tpb1 = tpb2

C
6662  FORMAT (L1,1X,A24,1X,I1,1X,4(2I3,1X))
6664  FORMAT (6E12.4)
6665  FORMAT (12(5L1,1X))
6666  FORMAT (12I6)
c slmod begin
6667  FORMAT (12I8)
c slmod end
66661 FORMAT (I3,1X,A6,1X,A4,A9,A3,2I3,3E12.4)
66662 FORMAT (2E12.4,10I6)
66664 FORMAT (I6,6X,5E12.4)
66665 FORMAT (2I6,4E12.4)
66666 FORMAT (I2,1X,A8,12(I3),1X,A10,1X,I2)
C
C   MODIFICATION OF INPUT DUE TO EITHER INCONSISTENCIES OR DUE
C   TO COUPLED NEUTRAL-PLASMA (OR NEUTRAL-NEUTRAL) CALCULATIONS
C   SOME FURTHER CONSTANTS ARE SET.      STATEM. NO. 2000 --> 4000
C
C
C   GEOMETRY, GRIDS
C
      IF (NLSLB) LEVGEO=1
      IF (NLCRC) LEVGEO=2
      IF (NLELL) LEVGEO=2
      IF (NLTRI) LEVGEO=2
      IF (NLPLG) LEVGEO=3
      IF (NLFEM) LEVGEO=4
      IF (NLTET) LEVGEO=5
      IF (NLGEN) LEVGEO=6
      IF (NLFEM) THEN
        NLPOL=.FALSE.
        NLTOR=.FALSE.
        IF (INDPRO(1).NE.3.AND.INDPRO(1).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(1),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' FINITE ELEMENT MESH '
          WRITE (iunout,*) 
     .      ' INDPRO(1) IS SET TO 3 <=> CONSTANT PROFILE '
          INDPRO(1)=3
        ENDIF
        IF (INDPRO(2).NE.3.AND.INDPRO(2).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(2),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' FINITE ELEMENT MESH '
          WRITE (iunout,*) 
     .      ' INDPRO(2) IS SET TO 3 <=> CONSTANT PROFILE '
          INDPRO(2)=3
        ENDIF
        IF (INDPRO(3).NE.3.AND.INDPRO(3).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(3),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' FINITE ELEMENT MESH '
          WRITE (iunout,*) 
     .      ' INDPRO(3) IS SET TO 3 <=> CONSTANT PROFILE '
          INDPRO(3)=3
        ENDIF
        IF (INDPRO(4).NE.3.AND.INDPRO(4).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(4),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' FINITE ELEMENT MESH '
          WRITE (iunout,*) 
     .      ' INDPRO(4) IS SET TO 3 <=> CONSTANT PROFILE '
          INDPRO(4)=3
        ENDIF
      ENDIF
      IF (NLTET) THEN
        IF (INDPRO(1).NE.3.AND.INDPRO(1).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(1),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' TETRAHEDRON MESH '
          WRITE (iunout,*) 
     .      ' INDPRO(1) IS SET TO 3 <=> CONSTANT PROFILE '
          INDPRO(1)=3
        ENDIF
        IF (INDPRO(2).NE.3.AND.INDPRO(2).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(2),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' TETRAHEDRON MESH '
          WRITE (iunout,*) 
     .      ' INDPRO(2) IS SET TO 3 <=> CONSTANT PROFILE '
          INDPRO(2)=3
        ENDIF
        IF (INDPRO(3).NE.3.AND.INDPRO(3).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(3),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' TETRAHEDRON MESH '
          WRITE (iunout,*) 
     .      ' INDPRO(3) IS SET TO 3 <=> CONSTANT PROFILE '
          INDPRO(3)=3
        ENDIF
        IF (INDPRO(4).NE.3.AND.INDPRO(4).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(4),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' TETRAHEDRON MESH '
          WRITE (iunout,*) 
     .      ' INDPRO(4) IS SET TO 3 <=> CONSTANT PROFILE '
          INDPRO(4)=3
        ENDIF
      ENDIF
C
      NP2ND=MIN0(N2ND,NP2ND)
      NT3RD=MIN0(N3RD,NT3RD)
      IF (NLTOR.AND.NLTRA) NTTRA=NT3RD
      NR1P2=0
      IF (NLPOL.OR.NLTOR) NR1P2=NR1ST
      NP2T3=0
      IF (NLTOR) NP2T3=NP2ND
      IF (.NOT.NLADD) NRADD=0
      IF (.NOT.NLMLT) NBMLT=1

      DO IN=1,NRAD
        NCLTAL(IN)=IN
      END DO
      NR1TAL = NR1ST
      NP2TAL = NP2ND
      NT3TAL = NT3RD
C
C  SOURCE PARAMETERS AND (REFLECTING) BOUNDARY CONDITIONS,
C  ON ADDITIONAL AND NON DEFAULT STANDARD SURFACES
C
      DO J=0,NLIMPS
C
        RINTEG(J)=RINTEG(1)
        EINTEG(J)=EINTEG(1)
        AINTEG(J)=AINTEG(1)
        DO ISPZ=1,NSPZ
          ISRS(ISPZ,J)=ISRS(1,J)
          ISRC(ISPZ,J)=ISRC(1,J)
          TRANSP(ISPZ,1,J)=TRANSP(1,1,J)
          TRANSP(ISPZ,2,J)=TRANSP(1,2,J)
          RECYCF(ISPZ,J)=RECYCF(1,J)
          RECYCT(ISPZ,J)=RECYCT(1,J)
          RECPRM(ISPZ,J)=RECPRM(1,J)
          EXPPL(ISPZ,J)=EXPPL(1,J)
          EXPEL(ISPZ,J)=EXPEL(1,J)
          EXPIL(ISPZ,J)=EXPIL(1,J)
          RECYCS(ISPZ,J)=RECYCS(1,J)
          RECYCC(ISPZ,J)=RECYCC(1,J)
          SPTPRM(ISPZ,J)=SPTPRM(1,J)
        end do
      end do

      DO WHILE (ASSOCIATED(REFLIST))
        NULLIFY(SURFCUR2)
        SURFCUR => SURFLIST
        DO WHILE (ASSOCIATED(SURFCUR))
          IF (SURFCUR%MODNAME == REFLIST%REFNAME) THEN
            NLJ = SURFCUR%NOSURF
            ILREF(NLJ) = REFLIST%JLREF
            ILSPT(NLJ) = REFLIST%JLSPT
            ISRS(:,NLJ) = REFLIST%JSRS
            ISRC(:,NLJ) = REFLIST%JSRC
            ZNML(NLJ) = REFLIST%ZNMLR
            EWALL(NLJ) = REFLIST%EWALLR
            EWBIN(NLJ) = REFLIST%EWBINR
            TRANSP(:,1,NLJ) = REFLIST%TRANSPR(:,1)
            TRANSP(:,2,NLJ) = REFLIST%TRANSPR(:,2)
            FSHEAT(NLJ) = REFLIST%FSHEATR
            RECYCF(:,NLJ) = REFLIST%RCYCFR
            RECYCT(:,NLJ) = REFLIST%RCYCTR
            RECPRM(:,NLJ) = REFLIST%RCPRMR
            EXPPL(:,NLJ) = REFLIST%EXPPLR
            EXPEL(:,NLJ) = REFLIST%EXPELR
            EXPIL(:,NLJ) = REFLIST%EXPILR
            RECYCS(:,NLJ) = REFLIST%RCYCSR
            RECYCC(:,NLJ) = REFLIST%RCYCCR
            SPTPRM(:,NLJ) = REFLIST%STPRMR
            IF (.NOT.ASSOCIATED(SURFCUR2)) THEN
              SURFLIST => SURFCUR%NEXT
              DEALLOCATE(SURFCUR)
              SURFCUR => SURFLIST
            ELSE
              SURFCUR2%NEXT => SURFCUR%NEXT
              DEALLOCATE(SURFCUR)
              SURFCUR => SURFCUR2%NEXT
            END IF
          ELSE
            SURFCUR2 => SURFCUR
            SURFCUR => SURFCUR%NEXT
          END IF
        END DO
        REFCUR => REFLIST
        REFLIST => REFLIST%NEXT
        DEALLOCATE (REFCUR%JSRS)
        DEALLOCATE (REFCUR%JSRC)
        DEALLOCATE (REFCUR%TRANSPR)
        DEALLOCATE (REFCUR%RCYCFR)
        DEALLOCATE (REFCUR%RCYCTR)
        DEALLOCATE (REFCUR%RCPRMR)
        DEALLOCATE (REFCUR%EXPPLR)
        DEALLOCATE (REFCUR%EXPELR)
        DEALLOCATE (REFCUR%EXPILR)
        DEALLOCATE (REFCUR%RCYCSR)
        DEALLOCATE (REFCUR%RCYCCR)
        DEALLOCATE (REFCUR%STPRMR)
        DEALLOCATE (REFCUR)
      ENDDO

      IF (ASSOCIATED(SURFLIST)) THEN
        WRITE (iunout,*) 
     .    ' SURFACE DATA HAVE NOT BEEN DEFINED FOR MODEL:'
        DO WHILE (ASSOCIATED(SURFLIST))
          WRITE (iunout,*) SURFLIST%MODNAME
          SURFCUR => SURFLIST
          SURFLIST => SURFLIST%NEXT
          DEALLOCATE(SURFCUR)
        END DO
        WRITE (iunout,*) ' EXECUTION IS STOPPED '
        CALL EXIT_OWN(1)
      END IF

      DO 2000 J=0,NLIMPS
        IF (ILCOL(J).LT.0) IGFIL(J)=1
        ILCOL(J)=MAX0(1,IABS(ILCOL(J)))
        IF (ILIIN(J).LE.0.OR.ILIIN(J).GE.3) ILSPT(J)=0
        IF (ILIIN(J).LE.0) TRANSP(:,1,J)=0.D0
        IF (ILIIN(J).LE.0) TRANSP(:,2,J)=0.D0
        ISPUT(1,J)=IDEZ(ILSPT(J),1,2)
        ISPUT(2,J)=IDEZ(ILSPT(J),2,2)
        IF (ILIIN(J).EQ.2) RECYCF(:,J)=0.
        IF (ILIIN(J).EQ.2) RECYCT(:,J)=0.
C
        SAVE=ZNML(J)
        ZNML(J)=DBLE(IDINT(SAVE/100.D0))
        ZNCL(J)=SAVE-100.*ZNML(J)
        DO 2001 ISPZ=1,NSPZ
          ISRF(ISPZ,J)=ISRF(ISPZ,1)
          ISRT(ISPZ,J)=ISRT(ISPZ,1)
2001    CONTINUE
2000  CONTINUE
C
      INMP1I=0
      INMP2I=0
      INMP3I=0
C
      DO 2019 ISTS=1,NSTSI
        NLJ=NLIM+ISTS
        ISWICH(1,NLJ)=IDEZ(ILSWCH(NLJ),1,6)
        IF (ISWICH(1,NLJ).EQ.1) ISWICH(1,NLJ)=-1
        IF (ISWICH(1,NLJ).EQ.2) ISWICH(1,NLJ)=1
        ISWICH(2,NLJ)=IDEZ(ILSWCH(NLJ),2,6)
        IF (ISWICH(2,NLJ).EQ.1) ISWICH(2,NLJ)=-1
        IF (ISWICH(2,NLJ).EQ.2) ISWICH(2,NLJ)=1
        ISWICH(3,NLJ)=IDEZ(ILSWCH(NLJ),3,6)
        IF (ISWICH(3,NLJ).EQ.1) ISWICH(3,NLJ)=-1
        IF (ISWICH(3,NLJ).EQ.2) ISWICH(3,NLJ)=1
        ISWICH(4,NLJ)=IDEZ(ILSWCH(NLJ),4,6)
        IF (ISWICH(4,NLJ).EQ.1) ISWICH(4,NLJ)=-1
        IF (ISWICH(4,NLJ).EQ.2) ISWICH(4,NLJ)=1
        ISWICH(5,NLJ)=IDEZ(ILSWCH(NLJ),5,6)
        IF (ISWICH(5,NLJ).EQ.1) ISWICH(5,NLJ)=-1
        IF (ISWICH(5,NLJ).EQ.2) ISWICH(5,NLJ)=1
        ISWICH(6,NLJ)=IDEZ(ILSWCH(NLJ),6,6)
        IF (ISWICH(6,NLJ).EQ.1) ISWICH(6,NLJ)=-1
        IF (ISWICH(6,NLJ).EQ.2) ISWICH(6,NLJ)=1
C
        IF (ISWICH(4,NLJ).NE.0.OR.ISWICH(5,NLJ).NE.0.OR.
     .      ISWICH(6,NLJ).NE.0) THEN
          ILBLCK(NLJ)=IDEZ(ILCELL(NLJ),4,4)
          I1000=1000*ILBLCK(NLJ)
          ILACLL(NLJ)=ILCELL(NLJ)-I1000
        ENDIF
C
        IF (ILSWCH(NLJ).NE.0.AND.ILIIN(NLJ).GT.0) THEN
          DO ISPZ=1,NSPZ
            IF (TRANSP(ISPZ,1,NLJ).NE.0..OR.TRANSP(ISPZ,2,NLJ).NE.0.)
     .      THEN
              WRITE (iunout,*) 
     .          'EXIT FROM TIMEA0: SURFACE NO. ISTS OPERATING'
              WRITE (iunout,*) 
     .          'A SWITCH BUT IS SOMETIMES TRANSPARENT AND '
              WRITE (iunout,*) 
     .          'SOMETIMES REFLECTING (SEMI-TRANSPARENCY OPTION)'
              WRITE (iunout,*) 
     .          'POSSIBLE FIXES: SEE MANUAL, CHAPTER 2, SECTION 6'
              WRITE (iunout,*) 'ISTS= ',ISTS
              CALL EXIT_OWN(1)
            ENDIF
          ENDDO
        ENDIF
C
C  SET NON DEFAULT STANDARD SURFACE IDENTIFIERS INMP...
C
C  RADIAL SURFACE
        DO 2014 IR=1,NR1ST
          IF (IR.EQ.INUMP(ISTS,1)) THEN
            INMP1I(IR,0,0)=ISTS
            DO J=IRPTA(ISTS,2),IRPTE(ISTS,2)-1
              INMP1I(IR,J,0)=ISTS
              DO K=IRPTA(ISTS,3),IRPTE(ISTS,3)-1
                INMP1I(IR,0,K)=ISTS
                INMP1I(IR,J,K)=ISTS
              END DO
            END DO
          ENDIF
2014    CONTINUE
C  POLOIDAL SURFACE
        DO 2016 JP=1,NP2ND
          IF (JP.EQ.INUMP(ISTS,2)) THEN
            INMP2I(0,JP,0)=ISTS
            DO I=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
              INMP2I(I,JP,0)=ISTS
              DO K=IRPTA(ISTS,3),IRPTE(ISTS,3)-1
                INMP2I(0,JP,K)=ISTS
                INMP2I(I,JP,K)=ISTS
              END DO
            END DO
          ENDIF
2016    CONTINUE
C  TOROIDAL SURFACE
        DO 2018 KT=1,NT3RD
          IF (KT.EQ.INUMP(ISTS,3)) THEN
            INMP3I(0,0,KT)=ISTS
            DO I=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
              INMP3I(I,0,KT)=ISTS
              DO J=IRPTA(ISTS,2),IRPTE(ISTS,2)-1
                INMP3I(0,J,KT)=ISTS
                INMP3I(I,J,KT)=ISTS
              END DO
            END DO
          ENDIF
2018    CONTINUE
C
2019  CONTINUE
C
      NLSYMT(0)=.TRUE.
      NLSYMP(0)=.TRUE.
      DO 2028 ISTRA=1,NSTRAI
        IF (INDSRC(ISTRA).EQ.6) GOTO 2028
        IF (.NOT.NLSRF(ISTRA))
     .  THMAX=MAX(0._DP,MIN(PIA,SORMAX(ISTRA)*DEGRAD))
        IF (NLSRF(ISTRA))
     .  THMAX=MAX(0._DP,MIN(PIHA,SORMAX(ISTRA)*DEGRAD))
        IF (NAMODS(ISTRA).EQ.1) THEN
          RP1=SORCOS(ISTRA)+1.
          SORCOS(ISTRA)=1./RP1
!pb          IF (ABS(COS(THMAX)).LE.EPS10) THEN
          IF (ABS(COS(THMAX)).LE.EPS5) THEN
            SORMAX(ISTRA)=1.
          ELSE
            SORMAX(ISTRA)=1.-COS(THMAX)**RP1
          ENDIF
        ELSEIF (NAMODS(ISTRA).EQ.2) THEN
          SORCOS(ISTRA)=SORCOS(ISTRA)*DEGRAD
          SORMAX(ISTRA)=THMAX
        ELSE
          WRITE (iunout,*) 'INPUT ERROR: ISTRA,NAMODS(ISTRA)='
          WRITE (iunout,*) '             ',ISTRA,NAMODS(ISTRA)
          CALL EXIT_OWN(1)
        ENDIF
        NLSYMT(0)=NLSYMT(0).AND.NLSYMT(ISTRA)
        NLSYMP(0)=NLSYMP(0).AND.NLSYMP(ISTRA)
2028  CONTINUE
C
C
C  SPECIES INDEX DISTRIBUTION OF PRIMARY SOURCE PARTICLES
C  OR FOR THERMAL PARTICLE REFLECTION MODEL
C
      NATMIM=NATMI-1
      NMOLIM=NMOLI-1
      NIONIM=NIONI-1
      NPLSIM=NPLSI-1
      NPHOTIM=NPHOTI-1
      SA=0.
      SI=0.
      SM=0.
      SPP=0.
      SPH=0.
      DO J=1,NIONI
        SI=SI+DIOD(J)
        DION(J)=SI
      END DO
      DO J=1,NMOLI
        SM=SM+DMLD(J)
        DMOL(J)=SM
      END DO
      DO J=1,NATMI
        SA=SA+DATD(J)
        DATM(J)=SA
      END DO
      DO J=1,NPLSI
        SPP=SPP+DPLD(J)
        DPLS(J)=SPP
      END DO
      DO J=1,NPHOTI
        SPH=SPH+DPHD(J)
        DPHOT(J)=SPH
      END DO
C
C  NORMALIZE DISTRIBUTION AND CUMULATIVE DISTRIBUTION
      DO IION=1,NIONI
        DIOD(IION)=DIOD(IION)/(SI+1.D-60)
        DION(IION)=DION(IION)/(SI+1.D-60)
      END DO
      DO IMOL=1,NMOLI
        DMLD(IMOL)=DMLD(IMOL)/(SM+1.D-60)
        DMOL(IMOL)=DMOL(IMOL)/(SM+1.D-60)
      END DO
      DO IATM=1,NATMI
        DATD(IATM)=DATD(IATM)/(SA+1.D-60)
        DATM(IATM)=DATM(IATM)/(SA+1.D-60)
      END DO
      DO IPLS=1,NPLSI
        DPLD(IPLS)=DPLD(IPLS)/(SPP+1.D-60)
        DPLS(IPLS)=DPLS(IPLS)/(SPP+1.D-60)
      END DO
      DO IPHOT=1,NPHOTI
        DPHD(IPHOT)=DPHD(IPHOT)/(SPH+1.D-60)
        DPHOT(IPHOT)=DPHOT(IPHOT)/(SPH+1.D-60)
      END DO
C
C
C  ATOMIC WEIGHT OF TEST IONS  =RMASSI(IION)
      DO IION=1,NIONI
        RMASSI(IION)=NMASSI(IION)*PMASSA
        RSQDVI(IION)=1._DP/SQRT(RMASSI(IION))*CVELAA
        CVRSSI(IION)=RMASSI(IION)*CVELI2
        ALMASI(IION)=LOG10(RMASSI(IION))
      END DO
C  ATOMIC WEIGHT OF ATOMS  =RMASSA(IATM)
      DO IATM=1,NATMI
        RMASSA(IATM)=NMASSA(IATM)*PMASSA
        RSQDVA(IATM)=1._DP/SQRT(RMASSA(IATM))*CVELAA
        CVRSSA(IATM)=RMASSA(IATM)*CVELI2
        ALMASA(IATM)=LOG10(RMASSA(IATM))
      END DO
C  ATOMIC WEIGHT OF MOLECULES
      DO IMOL=1,NMOLI
        RMASSM(IMOL)=NMASSM(IMOL)*PMASSA
        RSQDVM(IMOL)=1._DP/SQRT(RMASSM(IMOL))*CVELAA
        CVRSSM(IMOL)=RMASSM(IMOL)*CVELI2
        ALMASM(IMOL)=LOG10(RMASSM(IMOL))
      END DO
C  ATOMIC WEIGHT OF BULK IONS
      DO IPLS=1,NPLSI
        RMASSP(IPLS)=NMASSP(IPLS)*PMASSA
csw check photon in bulk
        if(rmassp(ipls) > 0._DP) then
          RSQDVP(IPLS)=1._DP/SQRT(RMASSP(IPLS))*CVELAA
          CVRSSP(IPLS)=RMASSP(IPLS)*CVELI2
          ALMASP(IPLS)=LOG10(RMASSP(IPLS))
        else
          rsqdvp(ipls)=0.
          cvrssp(ipls)=0.
          almasp(ipls)=0.
        endif
csw end check
      END DO
C
C
C  SET SOME ARRAYS TO SPEED UP COMPUTATIONS
C
C  1ST: SPECIES FLAGS:
        DO IPH=0,NPHOTP
        DO IA=0,NATMP
        DO IM=0,NMOLP
        DO II=0,NIONP
        DO IP=0,NPLSP
          ISPEZ(-1,IPH,IA,IM,II,IP)=-1
          ISPEZ(0,IPH,IA,IM,II,IP)=IPH
          ISPEZ(1,IPH,IA,IM,II,IP)=NPHOTI+IA
          ISPEZ(2,IPH,IA,IM,II,IP)=NPHOTI+NATMI+IM
          ISPEZ(3,IPH,IA,IM,II,IP)=NPHOTI+NATMI+NMOLI+II
          ISPEZ(4,IPH,IA,IM,II,IP)=NPHOTI+NATMI+NMOLI+NIONI+IP
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        DO IZ=1,NPHOTI
          ISPEZI(IZ,-1)=0
          ISPEZI(IZ,0)=IZ
          ISPEZI(IZ,1)=0
          ISPEZI(IZ,2)=0
          ISPEZI(IZ,3)=0
          ISPEZI(IZ,4)=0
        ENDDO
        DO IZ=NPHOTI+1,NPHOTI+NATMI
          ISPEZI(IZ,-1)=1
          ISPEZI(IZ,0)=0
          ISPEZI(IZ,1)=IZ-NPHOTI
          ISPEZI(IZ,2)=0
          ISPEZI(IZ,3)=0
          ISPEZI(IZ,4)=0
        ENDDO
        DO IZ=NPHOTI+NATMI+1,NPHOTI+NATMI+NMOLI
          ISPEZI(IZ,-1)=2
          ISPEZI(IZ,0)=0
          ISPEZI(IZ,1)=0
          ISPEZI(IZ,2)=IZ-NPHOTI-NATMI
          ISPEZI(IZ,3)=0
          ISPEZI(IZ,4)=0
        ENDDO
        DO IZ=NPHOTI+NATMI+NMOLI+1,NPHOTI+NATMI+NMOLI+NIONI
          ISPEZI(IZ,-1)=3
          ISPEZI(IZ,0)=0
          ISPEZI(IZ,1)=0
          ISPEZI(IZ,2)=0
          ISPEZI(IZ,3)=IZ-NPHOTI-NATMI-NMOLI
          ISPEZI(IZ,4)=0
        ENDDO
        DO IZ=NPHOTI+NATMI+NMOLI+NIONI+1,NPHOTI+NATMI+NMOLI+NIONI+NPLSI
          ISPEZI(IZ,-1)=4
          ISPEZI(IZ,0)=0
          ISPEZI(IZ,1)=0
          ISPEZI(IZ,2)=0
          ISPEZI(IZ,3)=0
          ISPEZI(IZ,4)=IZ-NPHOTI-NATMI-NMOLI-NIONI
        ENDDO

      IF (NPHOTI > 0) CALL PH_INIT(1)
      CALL SETAMD(0)
      CALL ALLOC_COMUSR(2)
      CALL ALLOC_CTEXT(2)

      CALL SETTXT

      IF (NADVI > 0) THEN
        TXTTAL(1:NADVI,NTALA) = TXTTLA(1:NADVI)
        TXTSPC(1:NADVI,NTALA) = TXTSCA(1:NADVI)
        TXTUNT(1:NADVI,NTALA) = TXTUTA(1:NADVI)
      END IF

      IF (NCLVI > 0) THEN
        TXTTAL(1:NCLVI,NTALC) = TXTTLC(1:NCLVI)
        TXTSPC(1:NCLVI,NTALC) = TXTSCC(1:NCLVI)
        TXTUNT(1:NCLVI,NTALC) = TXTUTC(1:NCLVI)
      END IF

      IF (NALVI > 0) THEN
        TXTTAL(1:NALVI,NTALR) = TXTTLR(1:NALVI)
        TXTSPC(1:NALVI,NTALR) = TXTSCR(1:NALVI)
        TXTUNT(1:NALVI,NTALR) = TXTUTR(1:NALVI)
      END IF

      IF (NSNVI > 0) THEN
        TXTTAL(1:NSNVI,NTALT) = TXTTLT(1:NSNVI)
        TXTSPC(1:NSNVI,NTALT) = TXTSCT(1:NSNVI)
        TXTUNT(1:NSNVI,NTALT) = TXTUTT(1:NSNVI)
      END IF

      DEALLOCATE (TXTTLA)
      DEALLOCATE (TXTSCA)
      DEALLOCATE (TXTUTA)
      DEALLOCATE (TXTTLC)
      DEALLOCATE (TXTSCC)
      DEALLOCATE (TXTUTC)
      DEALLOCATE (TXTTLR)
      DEALLOCATE (TXTSCR)
      DEALLOCATE (TXTUTR)
      IF (NPRNLI > 0) THEN
        DEALLOCATE (TXTTLT)
        DEALLOCATE (TXTSCT)
        DEALLOCATE (TXTUTT)
      END IF
C
C
C  ADDITIONAL INPUT FOR THIS RUN COMES FROM EITHER
C  ANOTHER CODE (DATA FILE) OR FROM AN EARLIER RUN OF EIRENE
C
3000  CONTINUE
C
C  INPUT BLOCK 14 BEGIN
C
C  READ DATA IN INTERFACING SUBROUTINE INFCOP  1400 -- 1499
C

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time vor block 14 ',tpb2-tpb1
!pb      tpb1 = tpb2

      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL MASAGE ('*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"   ')
      IF (NMODE.EQ.0) THEN
        WRITE (iunout,*) '        SUBR. INFCOP NOT CALLED. '
        READ (IUNIN,6666) NAINI,NCOPII,NCOPIE
        NCOPI=NCOPIE
        WRITE (iunout,*) '        NAINI, NCOPI = ',NAINI,NCOPI
        IF (NAINI.GT.NAIN) THEN
          CALL MASPRM('NAIN',4,NAIN,'NAINI',5,NAINI,IERROR)
          GOTO 1500
        ENDIF
        NCPVI=NCOPI*NPLSI
        IF (NCPVI.GT.NCPV) THEN
          CALL MASPRM('NCOP',4,NCOP,'NCOPI',5,NCOPI,IERROR)
          CALL EXIT_OWN(1)
        ENDIF
        CALL ALLOC_CCOUPL(2)
        DO 3020 J=1,NAINI
3021      READ (IUNIN,'(A72)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') GOTO 3021
          READ (ZEILE,6666) NAINS(J),NAINT(J)
          READ (IUNIN,'(A72)') TXTPLS(J,NTALN)
          READ (IUNIN,'(2A24)') TXTPSP(J,NTALN),TXTPUN(J,NTALN)
3020    CONTINUE
      ELSEIF (NMODE.NE.0) THEN
        NAINI=0
C  READ BLOCK 14 AND GEOMETRY FROM EXTERNAL DATABASE (FT30)
        CALL IF0COP
      ENDIF

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time nach if0cop ',tpb2-tpb1
!pb      tpb1 = tpb2

C
C  INPUT BLOCK 14 DONE
C
C
C  PLOTTING AND TEXT
C
      IF (NPLOTR.LE.0) NPLOTR=NR1ST
      IF (NPLOTP.LE.0) NPLOTP=NP2ND
      IF (NPLOTT.LE.0) NPLOTT=NT3RD
      NPLOTR=MIN0(NPLOTR,NR1ST)
      if (nltri.or.nltet) NPLOTR=MIN0(NPLOTR,NR1ST-1)
      NPLINR=MAX0(NPLINR,1)
      NPLDLR=MAX0(NPLDLR,1)
      NPLOTP=MIN0(NPLOTP,NP2ND)
      NPLINP=MAX0(NPLINP,1)
      NPLDLP=MAX0(NPLDLP,1)
      NPLOTT=MIN0(NPLOTT,NT3RD)
      NPLINT=MAX0(NPLINT,1)
      NPLDLT=MAX0(NPLDLT,1)
      PL1ST=PL1ST.AND.NR1ST.GT.1
      PL2ND=PL2ND.AND.NP2ND.GT.1
      PL3RD=PL3RD.AND.NT3RD.GT.1
      PLADD=PLADD.AND.NLIMI.GT.0
      NLPL2D=PL1ST.OR.PL2ND.OR.PL3RD.OR.PLADD
      NLPL3D=PL3A(1).OR.PL3A(2).OR.PL3A(3).OR.PL3A(4).OR.PL3A(5).OR.
     .       PL3S(1).OR.PL3S(2).OR.PL3S(3)
      PLHST=PLHST.AND.(NLPL2D.OR.NLPL3D)
      PLCHOR=PLCHOR.AND.NVOLPL.EQ.0
      TRCSUR=TRCSUR.AND.(NLIMI.GT.0.OR.NSTSI.GT.0)
      LPRADR=.FALSE.
      LPPOLR=.FALSE.
      LPTORR=.FALSE.

      IF ((NVOLPL > 0) .AND.
     .    (ANY(PLTL2D(1:NVOLPL)) .OR.
     .    (ANY(PLTL3D(1:NVOLPL)).AND.LEVGEO<=3))) CALL ALLOC_CGRPTL

      IF (NSIGVI > 0) THEN
        J=1
        DO WHILE (J <= NSIGVI)
          IF (LMISTALV(IIH(J))) THEN
            WRITE (iunout,*) 
     .        ' NO STATISTICS IS DONE FOR VOLUME TALLY NO. ',IIH(J)
            WRITE (iunout,*) ' BECAUSE TALLY HAS BEEN SWITCHED OFF '
            IGH(J)=IGH(NSIGVI)
            IIH(J)=IIH(NSIGVI)
            NSIGVI=NSIGVI-1
          ELSE
            J=J+1
          END IF
        END DO
      END IF

      IF (NSIGSI > 0) THEN
        J=1
        DO WHILE (J <= NSIGSI)
          IF (LMISTALS(IIHW(J))) THEN
            WRITE (iunout,*) 
     .        ' NO STATISTICS IS DONE FOR SURFACE TALLY NO. ',IIHW(J)
            WRITE (iunout,*) ' BECAUSE TALLY HAS BEEN SWITCHED OFF '
            IGHW(J)=IGHW(NSIGSI)
            IIHW(J)=IIHW(NSIGSI)
            NSIGSI=NSIGSI-1
          ELSE
            J=J+1
          END IF
        END DO
      END IF

      IF (NSIGCI > 0) THEN
        J=1
        DO WHILE (J <= NSIGCI)
          IF (LMISTALV(IIHC(1,J)) .OR. LMISTALV(IIHC(2,J))) THEN
            WRITE (iunout,*) 
     .        ' NO CORRELATION COEFFICIENT IS CALCULATED ',
     .        ' BETWEEN TALLIES ',IIHC(1,J),' AND ',IIHC(2,J)
            WRITE (iunout,*) ' BECAUSE TALLIES HAVE BEEN SWITCHED OFF '
            IGHC(1,J)=IGHC(1,NSIGCI)
            IIHC(1,J)=IIHC(1,NSIGCI)
            IGHC(2,J)=IGHC(2,NSIGCI)
            IIHC(2,J)=IIHC(2,NSIGCI)
            NSIGCI=NSIGCI-1
          ELSE
            J=J+1
          END IF
        END DO
      END IF

      IF (.NOT.LBGKV) NSIGI_BGK=0
      IF (.NOT.LCOPV) NSIGI_COP=0
C
C  NO MODIFICATION OF INPUT VARIABLES BEYOND THIS POINT
C  WITHOUT WARNING
C  EXCEPT IN SUBROUTINE MODUSR FOR THE NEXT ITERATION STEP
C         IN SUBROUTINE TMSUSR FOR THE NEXT TIME STEP
C
C
      CALL INIUSR
C
C  SET DERIVED INPUT PARAMETERS, GRIDS AND PROFILES
C
      NSURF=NR1ST*NP2ND*NT3RD*NBMLT
      NSTRD=NR1ST*NP2ND*NT3RD
      NBLCKS=NBMLT*NP2ND*NT3RD
      NSBOX=NSURF+NRADD
      NSBOX_TAL=NR1TAL*NP2TAL*NT3TAL*NBMLT+NRADD
      IF (NSBOX.GT.NRAD) THEN
        CALL MASPRM('NRAD',4,NRAD,'NSBOX',5,NSBOX,IERROR)
        CALL EXIT_OWN(1)
      ENDIF
      NSURFM=NSURF-1
      NR1STM=NR1ST-1
      NP2NDM=NP2ND-1
      NT3RDM=NT3RD-1
      NTTRAM=NTTRA-1
      NBMLTP=NBMLT+1
C
      NSIGI=NSIGVI+NSIGSI+NSIGCI
      NCPVI_STAT=0
      IF (NSIGI_COP > 0) NCPVI_STAT=NCPVI+NPLSI+2
      IF (NCPVI_STAT > NCPV_STAT) THEN
        CALL MASPRM('NCPVI_STAT',10,NCPVI_STAT,
     .              'NCPV_STAT',9,NCPV_STAT,IERROR)
        CALL EXIT_OWN(1)
      END IF
C
      CALL PAGE

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time vor grid ',tpb2-tpb1
!pb      tpb1 = tpb2

C
      IF (NFILEM.LE.1) THEN
C
C  SET GRIDS AND VOLUMES OF THE CELLS FOR THE STANDARD MESHES
C
C  SET RADIAL OR X GRID
        IF (NLRAD) CALL GRID (1)
C  SET POLOIDAL OR Y GRID
        IF (NLPOL) CALL GRID (2)
C  SET TOROIDAL OR Z GRID
        IF (NLTOR) CALL GRID (3)
C
        IF (INDPRO(12).LT.4) THEN
C  INITIALIZE SUBROUTINE VOLUME FOR LATER CALLS
C  TO BE WRITTEN:    CALL VOLUME(0)
C  SET VOLUMES, 1ST DIMENSION (R/X-GRID)
          IF (NLRAD) CALL VOLUME(1)
C  SET VOLUMES, 2ND DIMENSION (THETA/Y-GRID)
          IF (NLPOL) CALL VOLUME(2)
C  SET VOLUMES, 3RD DIMENSION (PHI/Z-GRID)
         IF (NLTOR) CALL VOLUME(3)
C  SET VOLUMES, IN ADDITIONAL CELL REGION
          IF (NLADD) CALL VOLUME(4)
        ELSEIF (INDPRO(12).EQ.4) THEN
C  READ VOLUMES FROM STREAM ISTREAM=VL0
          ISTREAM=VL0
          ITALI=NTALO
          CALL READTL(TXTPLS(1,ITALI),TXTPSP(1,ITALI),
     .                TXTPUN(1,ITALI),
     .                VOL,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .                3,ISTREAM)
C  TAKE VOLUMES FROM USER SUPPLIED ROUTINE
        ELSEIF (INDPRO(12).EQ.5) THEN
          CALL PROUSR (VOL,5+5*NPLS,0._DP,0._DP,0._DP,0._DP,
     .                              0._DP,0._DP,0._DP,NSBOX)
C  TAKE VOLUMES FROM EXTERNAL FILE
C  THIS CAN BE DONE ONLY AFTER CALL TO IF1COP
C  BECAUSE ONLY THEN WILL WORKARRAY BE FILLED FOR CELL VOLUMES
C       ELSEIF (INDPRO(12).EQ.6) THEN
C         CALL PROFR (VOL,5+5*NPLS,1,1,NSURF)
C         IF (NLADD) CALL VOLUME(4)
C       ELSEIF (INDPRO(12).EQ.7) THEN
C         CALL PROFR (VOL,5+5*NPLS,1,1,NSBOX)
        ENDIF
C
C  MULTIPLY GEOMETRICAL DATA, IF NBMLT.GT.1
C
        IF (NBMLT.GT.1) CALL MULTIG
C
C  SET CELL DIAMETER
C
        IF (LEVGEO == 5) THEN
!  Tetrahedra
c slmod begin - debug
c          WRITE(0,*) 'DEBUG: size celdia=',SIZE(celdia)
c          DO J = 1, NRAD
c            IF (VOL(J).GT.EPS10) THEN
c              CELDIA(J) = VOL(J)**(1._DP/3._DP)
c            ELSE
c              CELDIA(J) = 0._DP
c            ENDIF
c          ENDDO
c
          WHERE (VOL > EPS10)
            CELDIA = VOL**(1._DP/3._DP)
          ELSEWHERE
            CELDIA = 0._DP
          END WHERE
c slmod end
        ELSE
!  CYLINDRICAL OR TOROIDAL MESH
          WHERE (AREA > EPS10)
            CELDIA = SQRT(AREA/PIA)
          ELSEWHERE
            CELDIA = 0.D0
          END WHERE
        END IF
C
C   INCLUDE INFORMATION PROVIDED BY INPUT BLOCK 8: ADDITIONAL
C   DATA FOR SPECIFIC ZONES
C
        DO WHILE (ASSOCIATED(VOLLIST))
          VOL(VOLLIST%IN) = VOLLIST%VOL
          VOLCUR => VOLLIST
          VOLLIST => VOLLIST%NEXT
          DEALLOCATE(VOLCUR)
        ENDDO
        IF (NPHOTI > 0) CALL PH_INIT(2)
C
C   MODIFY SOME GEOMETRICAL DATA, USER SUPPLIED ROUTINE
C
        CALL GEOUSR
C
C  SET SOME DATA FOR ADDITIONAL SURFACES: INITIALIZE SUBR. TIMEA
C
        CALL TIMEA0
C
C   MODIFY THE BOUNDARIES OF SOME SURFACES TO AVOID ROUND OFF
C   ERRORS
C
        CALL SETFIT (TRCSUR)
C
C   SET THE COEFFICIENTS OF SOME SURFACES IDENTICAL TO THOSE
C   OF SOME OTHER, TO AVOID ROUND OFF ERRORS
C
        CALL SETEQ
C
        IF (LEVGEO.EQ.3.OR.LEVGEO.EQ.4) CALL WRMESH
C
C
        CALL INTVOL (VOL,1,1,NSBOX,VOLTOT,
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        WRITE (iunout,*) ' VOLTOT     ',VOLTOT
C
C  SET 'VISIBLE ADDITIONAL SURFACES' RANGES
Cc
        NLIMII=1
        NLIMIE=NLIMI
C
        NSOPT=MIN(NOPTIM,NSBOX)
        DO 8004 J=1,NSOPT
          IF (NLIMPB >= NLIMPS) THEN
            DO 8005 I=1,NLIMI
              LHELP(I) = IGJUM3(J,I)==0
8005        CONTINUE
            IIN=ILLZ(NLIMI,LHELP,1)+1
            IEN=NLIMI-ILLZ(NLIMI,LHELP,-1)
          ELSE
            IIN = 1 
            IEN = NLIMI 
! NO='1111....111'B ALL BITS SET TO 1
            NO=NOT(0)
! NUMBER OF INTEGERS USED TO STORE SURFACE INFORMATION
            IGO=NLIMI/NBITS 
            IF (MOD(NLIMI,NBITS) > 0) IGO = IGO + 1
! CHECK FOR FIRST ACTIVE SURFACE 
            DO I=1,IGO
              IF (IAND(IGJUM3(J,I),NO) /= NO) THEN
                IBEND = NBITS-1
                IF (I == IGO) IBEND = NLIMI-(I-1)*NBITS - 1
                DO IB=0,IBEND
                  IF (.NOT.BTEST(IGJUM3(J,I),IB)) THEN
                    ILA=IB
                    EXIT
                  END IF
                END DO
                IIN = (I-1)*NBITS+ILA+1
              END IF
            END DO
! CHECK FOR LAST ACTIVE SURFACE 
            DO I=IGO,1,-1
              IF (IAND(IGJUM3(J,I),NO) /= NO) THEN
                IBEND = NBITS-1
                IF (I == IGO) IBEND = NLIMI-(I-1)*NBITS - 1
                DO IB=IBEND,0,-1
                  IF (.NOT.BTEST(IGJUM3(J,I),IB)) THEN
                    ILA=IB
                    EXIT
                  END IF
                END DO
                IEN = (I-1)*NBITS+ILA+1
              END IF
            END DO
          END IF

          NLIMII(J)=IIN
          NLIMIE(J)=IEN
8004    CONTINUE
C
C  ALL GEOMETRICAL DATA (GRIDS, VOLUMES, SWITCHES) ARE DEFINED NOW
C
C   SAVE GEOMETRICAL DATA ON FILE FT12
C
        DO 8006 IRAD=1,NSBOX
          VOLG(IRAD)=VOL(IRAD)
8006    CONTINUE
C
        DO ILIMPS=1,NLMPGS
          AREAG(ILIMPS)=SAREA(ILIMPS)
        ENDDO
C
        IF (NFILEM.EQ.1) CALL WRGEOM(TRCFLE)
C
      ELSEIF (NFILEM.EQ.2) THEN
C
C   RESTORE GEOMETRICAL DATA FROM FILE FT12
C
        CALL RGEOM(TRCFLE)
C
        DO 8010 IRAD=1,NSBOX
          VOL(IRAD)=VOLG(IRAD)
8010    CONTINUE
C
        DO ILIMPS=1,NLMPGS
          SAREA(ILIMPS)=AREAG(ILIMPS)
        ENDDO
C
C
      ENDIF
C
      DO 8011 IS=1,NLIMPS
        IF (ILACLL(IS).NE.0.AND..NOT.NLADD) THEN
          WRITE (iunout,*) 'ADDITIONAL CELL SWITCHES DEFINED, BUT NO '
          WRITE (iunout,*) 'ADDITIONAL CELLS DEFINED'
          ISS=IS
          IF (ISS.GT.NLIM) ISS=-(ISS-NLIM)
          WRITE (iunout,*) 'SURFACE NO. IS= ',ISS
          WRITE (iunout,*) 'CHECK INPUT BLOCK 2E'
          CALL EXIT_OWN(1)
        ELSEIF (ILBLCK(IS).NE.0.AND..NOT.(NLMLT.OR.NLADD)) THEN
          WRITE (iunout,*) 'BLOCK SWITCHES DEFINED, BUT NEITHER BLOCKS'
          WRITE (iunout,*) 'NOR ADDITIONAL CELLS DEFINED'
          ISS=IS
          IF (ISS.GT.NLIM) ISS=-(ISS-NLIM)
          WRITE (iunout,*) 'SURFACE NO. IS= ',ISS
          WRITE (iunout,*) 'CHECK INPUT BLOCKS 2D AND 2E'
          CALL EXIT_OWN(1)
        ENDIF
8011  CONTINUE
C
4000  CONTINUE
C

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time vor plasma definition ',tpb2-tpb1
!pb      tpb1 = tpb2

      IF (ANY(INDPRO(1:6) == 6)) CALL ALLOC_BCKGRND
      IF (NMODE.NE.0.AND.IITER.LE.1) THEN
C  READ PLASMA BACKGROUND FROM EXTERNAL DATABASE (FT31) (NOT NLPLAS)
C  OR FROM COMMON BRAEIR (NLPLAS)
        CALL IF1COP
      ENDIF

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time fuer if1cop ',tpb2-tpb1
!pb      tpb1 = tpb2
C
      IF (NSTEP > 0) CALL ALLOC_CSTEP
C
!PB VOLTAL INITIALIZED WITH EPS60 TO AVOID ZERODIVISIONS (E.G. IN CUT CELLS)
      VOLTAL = EPS60
      DO IN=1,NSBOX
        INT = NCLTAL(IN)
        IF (INT > 0) VOLTAL(INT) = VOLTAL(INT) + VOL(IN)
      END DO
      CALL INTVOL (VOLTAL,1,1,NSBOX_TAL,VOLTOT_TAL,
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      WRITE (iunout,*) ' VOLTOT_TAL ',VOLTOT_TAL

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time fuer intvol(voltal) ',tpb2-tpb1
!pb      tpb1 = tpb2
C
      IF ((NFILEL.LE.1).OR.(NFILEL == 6)) THEN
C
C  SET PLASMA PARAMETERS AND SOURCE PARAMETERS
C
        CALL PLASMA

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time fuer plasma ',tpb2-tpb1
!pb      tpb1 = tpb2
C
C  MULTIPLY PLASMA PARAMETERS, IF NBLCKS.GT.1
C
        IF (NBLCKS.GT.1) CALL MULTIP
C
C   INCLUDE INFORMATION PROVIDED BY INPUT BLOCK 8: ADDITIONAL
C   DATA FOR SPECIFIC ZONES
C
        IF (NZADD.GT.0) THEN
          DO WHILE (ASSOCIATED(TEMPLIST))
            TEIN(TEMPLIST%IN) = TEMPLIST%TE
            IPLSTI = MPLSTI(TEMPLIST%IDION)
            TIIN(IPLSTI,TEMPLIST%IN) = TEMPLIST%TI
            TEMPCUR => TEMPLIST
            TEMPLIST => TEMPLIST%NEXT
            DEALLOCATE(TEMPCUR)
          ENDDO

          DO WHILE (ASSOCIATED(DENLIST))
            DIIN(DENLIST%IDION,DENLIST%IN) = DENLIST%DI
            DENCUR => DENLIST
            DENLIST => DENLIST%NEXT
            DEALLOCATE(DENCUR)
          ENDDO

          DO WHILE (ASSOCIATED(VELLIST))
            IPLS = VELLIST%IDION
            IPLSTI = MPLSTI(IPLS)
            IPLSV = MPLSV(IPLS)
            J = VELLIST%IN
            VXIN(IPLSV,J) = VELLIST%VX
            VYIN(IPLSV,J) = VELLIST%VY
            VZIN(IPLSV,J) = VELLIST%VZ
            IF (VELLIST%IZ .EQ. 1) THEN
              VXIN(IPLSV,J)=CVEL2A*VXIN(IPLSV,J)*
     .                      SQRT((TEIN(J)+TIIN(IPLSTI,J))/RMASSP(IPLS))
              VYIN(IPLSV,J)=CVEL2A*VYIN(IPLSV,J)*
     .                      SQRT((TEIN(J)+TIIN(IPLSTI,J))/RMASSP(IPLS))
              VZIN(IPLSV,J)=CVEL2A*VZIN(IPLSV,J)*
     .                      SQRT((TEIN(J)+TIIN(IPLSTI,J))/RMASSP(IPLS))
            ENDIF
            VELCUR => VELLIST
            VELLIST => VELLIST%NEXT
            DEALLOCATE(VELCUR)
          ENDDO
        ENDIF
C
C  MODIFY SOME PLASMA DATA, USER SUPPLIED ROUTINE
C
        CALL PLAUSR

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time fuer plausr ',tpb2-tpb1
!pb      tpb1 = tpb2

C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE INPUT PROFILES
C
        CALL PLASMA_DERIV(0)

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time fuer plasma_deriv ',tpb2-tpb1
!pb      tpb1 = tpb2
C
C  SET ATOMIC DATA TABLES
C
        CALL SETAMD(1)

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time fuer setamd ',tpb2-tpb1
!pb      tpb1 = tpb2
C
        IF (NFILEL.EQ.1) CALL WRPLAM(TRCFLE,0)
        IF (NFILEL.EQ.6) CALL WRPLAM_XDR(TRCFLE,0)

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time fuer wrplam ',tpb2-tpb1
!pb      tpb1 = tpb2
C
      ELSEIF (NFILEL.EQ.2.OR.NFILEL.EQ.3.OR.NFILEL.EQ.4.OR.
     .        NFILEL.EQ.7.OR.NFILEL.EQ.8.OR.NFILEL.EQ.9) THEN
C
C  READ PLASMA DATA, ATOMIC DATA, SOURCE DATA FROM FT13
C
        IF ((NFILEL == 2) .OR. (NFILEL == 3)) THEN
          CALL RPLAM(TRCFLE,0)
        ELSEIF ((NFILEL == 7) .OR. (NFILEL == 8)) THEN
          CALL RPLAM_XDR(TRCFLE,0)
        ELSEIF (NFILEL == 4) THEN
          CALL RPLAM(TRCFLE,NFILEL)
        ELSEIF (NFILEL == 9) THEN
          CALL RPLAM_XDR(TRCFLE,NFILEL)
        END IF
        CALL XSECTPH
C
      ENDIF

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time nach plasma definition ',tpb2-tpb1
!pb      tpb1 = tpb2

C
C  SETUP TABLE OF CONTRIBUTIONS OF MONTE-CARLO PARTICLES TO BACKGROUND SPECIES
C
      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)
      
      DO IPLS = 1, NPLSI
        IF (LEN_TRIM(CDENMODEL(IPLS)) > 0) THEN
          DO IRE = 1, TDMPAR(IPLS)%TDM%NRE
            ITYP = TDMPAR(IPLS)%TDM%ITP(IRE)
            ISPZ = IADTYP(ITYP) + TDMPAR(IPLS)%TDM%ISP(IRE)
            ISPZ_BACK(ISPZ,IPLS) = 1
          END DO
        END IF      
      END DO 

      CALL LEER(2)
      WRITE (IUNOUT,*) ' LIST OF CONTRIBUTIONS TO BACKGROUND SPECIES '
      DO ISPZ = 1, NSPZ
        DO IPLS = 1, NPLSI
          IF (ISPZ_BACK(ISPZ,IPLS) > 0)
     .      WRITE (IUNOUT,*) TEXTS(ISPZ), ' CONTRIBUTES TO ',
     .                       TEXTS(NSPAMI+IPLS)
        END DO
      END DO
C
C  AT THIS POINT THE BACKGROUND MEDIUM DATA ARE ALL SET.
C
C  COMPUTE SOURCE DATA (OVERRULE SOME OF INPUT BLOCK 7)
      IF ((NMODE.NE.0.AND.IITER.LE.1) .OR. (IITER > NITER))  THEN
        DO ISTRA=1,NSTRAI
          IF (INDSRC(ISTRA).GE.0) CALL IF2COP(ISTRA)
        ENDDO
      ENDIF
C
C
      IF (TRCSUR) THEN
C
        CALL LEER(2)
        WRITE (iunout,*) 'COEFFICIENTS FOR ADDITIONAL SURFACES'
        WRITE (iunout,*) 
     .    'THIS IS AFTER GEOUSR, SETEQ AND SETFIT ARE CALLED '
        DO 7701 J=1,NLIMI
          CALL LEER(2)
          WRITE (iunout,*) TXTSFL(J)
          CALL LEER(1)
          IF (IGJUM0(J) == 1) THEN
            WRITE (iunout,*) 'THIS SURFACE IS NOT DEFINED'
          ELSE
            WRITE (iunout,*) 'A0       ',A0LM(J)
            WRITE (iunout,*) 'A1,A2,A3 ',A1LM(J),A2LM(J),A3LM(J)
            WRITE (iunout,*) 'A4,A5,A6 ',A4LM(J),A5LM(J),A6LM(J)
            WRITE (iunout,*) 'A7,A8,A9 ',A7LM(J),A8LM(J),A9LM(J)
            WRITE (iunout,*) 'P:', P1(1,j),P1(2,j),P1(3,j),P2(1,j),
     .                             P2(2,j),P2(3,j)
            WRITE (iunout,*) 'JUMLIM ',JUMLIM(J)
            WRITE (iunout,*) 'ISWICH(1),ISWICH(2),ISWICH(3),ISWICH(4),',
     .                  'ISWICH(5),ISWICH(6)'
            WRITE (iunout,*)  ISWICH(1,J),ISWICH(2,J),ISWICH(3,J),
     .                   ISWICH(4,J),ISWICH(5,J),ISWICH(6,J)
          ENDIF
7701    CONTINUE
        CALL LEER(2)
        CALL MASIR2('IGJUM0 ',IGJUM0,1,1,1,1,NLIMPS)
        CALL LEER(1)
        IF (NLIMPB >= NLIMPS) THEN
          CALL MASIR2('IGJUM1 ',IGJUM1,0,NLIMPS,1,NLIMPS,NLIMPS)
        ELSE
          CALL MASBR2('IGJUM1 ',IGJUM1,0,NLIMPS,1,NLIMPS,NLIMPS,NBITS)
        END IF
        CALL LEER(1)
        IF (NLIMPB >= NLIMPS) THEN
          CALL MASIR2('IGJUM2 ',IGJUM2,0,NLIMPS,1,NLIMPS,NLIMPS)
        ELSE
          CALL MASBR2('IGJUM2 ',IGJUM2,0,NLIMPS,1,NLIMPS,NLIMPS,NBITS)
        END IF
        CALL LEER(1)
        NSOPT=MIN(NSBOX,NOPTIM)
        IF (NLIMPB >= NLIMPS) THEN
          CALL MASIR2('IGJUM3 ',IGJUM3,0,NOPTIM,1,NSOPT,NLIMPS)
        ELSE
          CALL MASBR2('IGJUM3 ',IGJUM3,0,NOPTIM,1,NSOPT,NLIMPS,NBITS)
        END IF
        CALL LEER(1)
        DO 7702 J=1,NSOPT
          CALL MASJ3 ('J,NLIMII,NLIMIE          ',J,NLIMII(J),NLIMIE(J))
7702    CONTINUE
C
      ENDIF

      CALL DEALLOC_BCKGRND
C
      IF (NPHOTI > 0) CALL PH_INIT(3)

!  allocate and initialize storage for trajectories 
      IF (.NOT.ALLOCATED(TRAJ)) THEN
        ALLOCATE (TRAJ(NCHORI+NTRJ))

        DO ITRJ = 1, NCHORI+NTRJ
          ALLOCATE(TRAJ(ITRJ)%TRJ)
          TRAJ(ITRJ)%TRJ%NCOU_CELL = 0
          NULLIFY(TRAJ(ITRJ)%TRJ%CELLS)
        END DO
      END IF

      IF (NCHORI > 0) THEN
        IF (ANY(NLSTCHR)) CALL SETUP_CHORD_SPECTRA
      END IF

!  determine number of background spectra

      NBACK_SPEC = 0

      DO J = 1, NADSPC
!  spectrum in geometrical cell 
        IF (ESTIML(J)%PSPC%ISRFCLL == 2) THEN
          ISPZ=IADTYP(ESTIML(J)%PSPC%IPRTYP) + ESTIML(J)%PSPC%IPRSP 
          NBACK_SPEC = NBACK_SPEC + COUNT(ISPZ_BACK(ISPZ,:)>0)
          LSPCCLL(ESTIML(J)%PSPC%ISPCSRF) = .TRUE.
        END IF
      END DO 

      IF (.NOT.ALLOCATED(BACK_SPEC) .AND. (NBACK_SPEC > 0)) 
     .   ALLOCATE(BACK_SPEC(NBACK_SPEC))

!pb      TPB2=SECOND_OWN()     
!pb      write (iunout,*) ' cpu-time am ende von input ',tpb2-tpb1
!pb      tpb1 = tpb2

C
      RETURN
C
C  ERROR EXITS
C
990   CONTINUE
      WRITE (iunout,*) 'TALLY NUMBER FOR PRINTOUT OR PLOT OF VOLUME'
      WRITE (iunout,*) 'AVERAGED TALLIES OUT OF RANGE'
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'TALLY NUMBER FOR PRINTOUT OF SURFACE AVERAGED '
      WRITE (iunout,*) 'TALLIES OUT OF RANGE'
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'FINITE ELEMENT OPTION USED, BUT GRID INDICATOR '
      WRITE (iunout,*) 'LESS THAN 6. '
      CALL EXIT_OWN(1)
994   CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT: NRPLG.NE.NP2ND, BUT NLPOL=TRUE'
      WRITE (iunout,*) 'NRPLG,NP2ND ',NRPLG,NP2ND
      CALL EXIT_OWN(1)
995   CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT: ',
     .            'WRONG UNIT NUMBER FOR OUTPUT OF TALLY SPECIFIED'
      WRITE (iunout,*) 'NTLV,NTLVF',NTLV,NTLVF
      CALL EXIT_OWN(1)
998   CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT BLOCK FOR ADDITIONAL DATA FOR  '
      WRITE (iunout,*) 'SPECIFIC ZONES FOUND AT ZONE NO. ',I
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: mkcens.f

C
      SUBROUTINE MKCENS

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CINIT
      USE COMPRT
      USE COMNNL
      USE COMSOU
      USE CTEXT

      IMPLICIT NONE

      REAL(DP) :: DTIMVO
      INTEGER :: I, IFIRST, IERROR

      SAVE
      DATA IFIRST/0/
C
C  SET DEFAULTS FOR SOURCE DUE TO INITIAL CONDITION, VALID ONLY FOR
C  FIRST TIMESTEP. MODIFIED FOR LATER TIMESTEPS IN SUBR. TMSTEP
C
C  DEFINE ONE MORE STRATUM
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        NSTRAI=NSTRAI+1
C  CHECK STORAGE
        IF (NSTRAI.GT.NSTRA) THEN
          CALL MASPRM('NSTRA',5,NSTRA,'NSTRAI',6,NSTRAI,IERROR)
          CALL EXIT_OWN(1)
        ENDIF
C
        TXTSOU(NSTRAI)='SOURCE DUE TO INITIAL CONDITION          '
C  SOURCE DISTRIBUTION SAMPLED FROM CENSUS ARRAYS RPARTC,IPARTC
        NLCNS(NSTRAI)=.TRUE.
        NLPNT(NSTRAI)=.FALSE.
        NLLNE(NSTRAI)=.FALSE.
        NLSRF(NSTRAI)=.FALSE.
        NLVOL(NSTRAI)=.FALSE.
C  DO NOT CALL IF2COP(NSTRAI)
        INDSRC(NSTRAI)=-1
C
        NLAVRP(NSTRAI)=.FALSE.
        NLAVRT(NSTRAI)=.FALSE.
        NLSYMP(NSTRAI)=.FALSE.
        NLSYMT(NSTRAI)=.FALSE.
        NPTS(NSTRAI)=0
        NINITL(NSTRAI)=2000*NINITL(NSTRAI-1)+1
        NEMODS(NSTRAI)=1
        NAMODS(NSTRAI)=1
        FLUX(NSTRAI)=0.
        NLATM(NSTRAI)=.FALSE.
        NLMOL(NSTRAI)=.FALSE.
        NLION(NSTRAI)=.FALSE.
        NLPLS(NSTRAI)=.FALSE.
        NSPEZ(NSTRAI)=0
        NSRFSI(NSTRAI)=0
C
        SORENI(NSTRAI)=0.
        SORENE(NSTRAI)=0.
        SORVDX(NSTRAI)=0.
        SORVDY(NSTRAI)=0.
        SORVDZ(NSTRAI)=0.
        SORCOS(NSTRAI)=0.
        SORMAX(NSTRAI)=0.
        SORCTX(NSTRAI)=0.
        SORCTY(NSTRAI)=0.
        SORCTZ(NSTRAI)=0.
C
      ENDIF
C
C  READ INITIAL POPULATION FROM PREVIOUS RUN, OVERWRITE DEFAULTS
C
      IF (NFILEJ.EQ.2.OR.NFILEJ.EQ.3) THEN

C  NEW TIMESTEP
        IF (DTIMVN.LE.0.D0) THEN
          DTIMVN=DTIMV
C       ELSE
C         DTIMVN=DTIMVN
        ENDIF
C
C  READ CENSUS ARRAY FROM OLD TIMESTEP
        CALL RSNAP
        DTIMVO=DTIMV
C
        WRITE (iunout,*) 'INITIAL POPULATION FOR FIRST TIMESTEP'
        WRITE (iunout,*) 'READ FROM FILE FT 15 '
        WRITE (iunout,*) 'PARTICLES AND FLUX STORED FOR INITIAL '
        WRITE (iunout,*) 'DISTRIBUTION IN PREVIOUS RUN '
        CALL MASJ1('IPRNL   ',IPRNL)
        CALL MASR1('FLUX    ',FLUX(NSTRAI))
C
        IF (DTIMVN.NE.DTIMVO) THEN
          FLUX(NSTRAI)=FLUX(NSTRAI)*DTIMVO/DTIMVN
C
          WRITE (iunout,*) 'FLUX IS RESCALED BY DTIMV_OLD/DTIMV_NEW '
          CALL MASR1('FLUX    ',FLUX(NSTRAI))
          CALL LEER(1)
        ENDIF
C
        DTIMV=DTIMVN
C
        IF (NPTST.LE.0) THEN
          NPTS(NSTRAI)=IPRNL
        ELSEIF (NPTST.GT.0) THEN
          NPTS(NSTRAI)=NPTST
        ENDIF
C
        IF (NPTS(NSTRAI).GT.0.AND.FLUX(NSTRAI).GT.0) THEN
          NSRFSI(NSTRAI)=1
          SORWGT(1,NSTRAI)=1.D0
        ENDIF
C
        CALL LEER(2)
        DO I=1,IPRNL
          RPARTC(I,10)=TIME0
        ENDDO
        WRITE (iunout,*) 'PARTICLE CLOCK RESET TO TIME0'
        WRITE (iunout,*) 'FIRST TIMESTEP RUNS FROM TIM1 TO TIM2:  '
          CALL MASR2('TIM1, TIM2      ',TIME0,TIME0+DTIMV)
          CALL LEER(2)
C
      ENDIF
      RETURN
      END
C ===== SOURCE: multi.f
C
C  may 05:  "no multip on averaging cells" corrected for 3D grids
C
      SUBROUTINE MULTI

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CINIT
      USE CGRID
      USE CGEOM

      IMPLICIT NONE

      INTEGER :: I, J, K
C
C  GEOMETRY DATA
C
      ENTRY MULTIG
C
C  ZONE VOLUMES, KNOWN IN ZONE 1 TO NSTRD
      DO 130 J=2,NBMLT
        DO 120 I=1,NSTRD
          VOL(I+(J-1)*NSTRD)=VOL(I)*VOLCOR(J)
120     CONTINUE
130   CONTINUE
      DO 140 I=1,NSTRD
        VOL(I)=VOL(I)*VOLCOR(1)
140   CONTINUE
C
      RETURN
C
C  PLASMA DATA
C
      ENTRY MULTIP
C
C  INDPRO.LE.4: ONLY RADIAL PLASMA PROFILES GIVEN
C  RADIAL PLASMA PROFILES, KNOWN IN ZONES 1 TO NR1ST
C  NBLCKS=NP2ND*NT3RD*NBMLT
C
      DO 210 J=2,NBLCKS
C  RADIAL BLOCK NO J
C  IS THIS A SPACE FOR AVERAGING: THEN DO NOT COPY
        IF (NP2ND.GT.1.AND.MOD(J,NP2ND).EQ.0) GOTO 210
        IF (NT3RD.GT.1.AND.NP2ND.LE.1.AND.MOD(J,NT3RD).EQ.0) GOTO 210
        IF (NT3RD.GT.1.AND.NP2ND.GT.1.AND.J.GT.NP2ND*(NT3RD-1)) GOTO 210
        IF (INDPRO(1).LE.4) THEN
          DO 201 I=1,NR1ST
            TEIN(I+(J-1)*NR1ST)=TEIN(I)
201       CONTINUE
        ENDIF
        IF (INDPRO(2).LE.4) THEN
          DO 202 K=1,NPLSTI
            DO 202 I=1,NR1ST
              TIIN(K,I+(J-1)*NR1ST)=TIIN(K,I)
202       CONTINUE
        ENDIF
        IF (INDPRO(3).LE.4) THEN
          DO 204 K=1,NPLSI
            DO 204 I=1,NR1ST
              DIIN(K,I+(J-1)*NR1ST)=DIIN(K,I)
204       CONTINUE
        ENDIF
        IF (INDPRO(4).LE.4) THEN
          DO 205 K=1,NPLSV
            DO 205 I=1,NR1ST
              VXIN(K,I+(J-1)*NR1ST)=VXIN(K,I)
              VYIN(K,I+(J-1)*NR1ST)=VYIN(K,I)
              VZIN(K,I+(J-1)*NR1ST)=VZIN(K,I)
205       CONTINUE
        ENDIF
        IF (INDPRO(6).LE.4) THEN
          DO 207 K=1,NAINI
          DO 207 I=1,NR1ST
            ADIN(K,I+(J-1)*NR1ST)=ADIN(K,I)
207       CONTINUE
        ENDIF
210   CONTINUE
C
C  INDPRO.GT.4: ONLY NSTRD=NR1ST*NP2ND*NT3RD PLASMA DATA GIVEN
      DO 310 J=2,NBMLT
        IF (INDPRO(1).GT.4) THEN
          DO 301 I=1,NSTRD
            TEIN(I+(J-1)*NSTRD)=TEIN(I)
301       CONTINUE
        ENDIF
        IF (INDPRO(2).GT.4) THEN
          DO 302 K=1,NPLSTI
            DO 302 I=1,NSTRD
              TIIN(K,I+(J-1)*NSTRD)=TIIN(K,I)
302       CONTINUE
        ENDIF
        IF (INDPRO(3).GT.4) THEN
          DO 304 K=1,NPLSI
            DO 304 I=1,NSTRD
              DIIN(K,I+(J-1)*NSTRD)=DIIN(K,I)
304       CONTINUE
        ENDIF
        IF (INDPRO(4).GT.4) THEN
          DO 305 K=1,NPLSV
            DO 305 I=1,NSTRD
              VXIN(K,I+(J-1)*NSTRD)=VXIN(K,I)
              VYIN(K,I+(J-1)*NSTRD)=VYIN(K,I)
              VZIN(K,I+(J-1)*NSTRD)=VZIN(K,I)
305       CONTINUE
        ENDIF
        IF (INDPRO(5).GT.4) THEN
          DO 306 I=1,NSTRD
            BXIN(I+(J-1)*NSTRD)=BXIN(I)
            BYIN(I+(J-1)*NSTRD)=BYIN(I)
            BZIN(I+(J-1)*NSTRD)=BZIN(I)
306       CONTINUE
        ENDIF
        IF (INDPRO(6).LE.4) THEN
          DO 307 K=1,NAINI
          DO 307 I=1,NSTRD
            ADIN(K,I+(J-1)*NSTRD)=ADIN(K,I)
307       CONTINUE
        ENDIF
310   CONTINUE
C
      RETURN
      END
C ===== SOURCE: nanalg.f
C
C   ****************************
C   *  NONANALOG METHODS, A.I. *
C   ****************************
C
C      SUBROUTINE NANALG
C
      SUBROUTINE NANALG
C
C   SET UP SPLITTING SURFACES
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE CLOGAU
      USE CGRID
      USE CTRCEI
      USE COMPRT
      USE COMSPL

      IMPLICIT NONE

      REAL(DP) :: RSPLIT(N1ST)
      REAL(DP) :: XMAXR1, ZR1, ZR2
      INTEGER :: IR, IP, IT, IRT, IRP, IA, IRA, IRD, J, MMXRAD, MAXR1,
     .           JR, JS

C---------------------------------------------------------------------
      NODES=0
      RSPLIT=0.D0
      NLSPLT(1:N1ST+N2ND+N3RD+NLIM)=.FALSE.
C
C
100   CONTINUE
C  SET RADIAL SPLITTING SURFACES
      IF(NR1ST.LE.2) GOTO 200
C
C  DEFAULT MODEL: 101--119
C
C  HARD WIRED OPTIONS FOR RADIAL SPLITTING SURFACES
C  USES: MAXRAD, SPLPAR, AND RHOSRF-GRID
C  SETS: RSPLIT(ISURF) : SPLITTING SURFACE RADIA
C  FINDS: JS: RADIAL SURFACE NUMBER TO BE USED AS SPLITTING SURFACE
C
      IF (LEVGEO.LE.3.AND.MAXRAD.LT.0) THEN
        MMXRAD=ABS(MAXRAD)
        MAXR1=MMXRAD+1
        XMAXR1=MMXRAD+1
        DO 101 J=2,MAXR1
          RSPLIT(J-1)=RHOSRF(1)+(J-1)/XMAXR1*(RHOSRF(NR1ST)-RHOSRF(1))
101     CONTINUE
C
C
C
        DO 110 JS=2,MIN(NR1ST,N1ST+N2ND+N3RD+NLIM)
          ZR1=RHOSRF(JS-1)
          ZR2=RHOSRF(JS)
          DO 104 JR=1,MMXRAD
            IF(RSPLIT(JR).GE.ZR2.OR.RSPLIT(JR).LT.ZR1) GOTO 104
C  RADIAL SURFACE JS SHOULD BE USED AS SPLITTING SURFACE
            RNUMB(JS)=SPLPAR
            NLSPLT(JS)=.TRUE.
104       CONTINUE
110     CONTINUE
        GOTO 200
      ENDIF
C
C  NON DEFAULT MODEL: 121--199
C
      DO 121 IRD=1,MAXRAD
        DO IR=1,MIN(NR1ST,N1ST+N2ND+N3RD+NLIM)
          IF (IR.EQ.NSSPL(IRD)) THEN
            RNUMB(IR)=PRMSPL(IRD)
            NLSPLT(IR)=.TRUE.
          ENDIF
        ENDDO
121   CONTINUE
C
200   CONTINUE
C
C  SET POLOIDAL SPLITTING SURFACES
      IF(NP2ND.LE.2) GOTO 300
C
C  NON DEFAULT MODEL: 221--299
C
      DO 221 IRP=1,MAXPOL
        DO IP=1,NP2ND
          IF (IP.EQ.NSSPL(N1ST+IRP)) THEN
            RNUMB(N1ST+IP)=PRMSPL(N1ST+IRP)
            NLSPLT(N1ST+IP)=.TRUE.
          ENDIF
        ENDDO
221   CONTINUE
C
300   CONTINUE
C
C  SET TOROIDAL SPLITTING SURFACES
      IF(NT3RD.LE.2) GOTO 400
C
C  NON DEFAULT MODEL: 321--399
C
      DO 321 IRT=1,MAXTOR
        DO IT=1,NT3RD
          IF (IT.EQ.NSSPL(N1ST+N2ND+IRT)) THEN
            RNUMB(N1ST+N2ND+IT)=PRMSPL(N1ST+N2ND+IRT)
            NLSPLT(N1ST+N2ND+IT)=.TRUE.
          ENDIF
        ENDDO
321   CONTINUE
C
400   CONTINUE
C  SET ADDITIONAL SPLITTING SURFACES
      IF(NLIMI.LE.0) GOTO 500
C
C  NON DEFAULT MODEL: 421--499
C
      DO 421 IRA=1,MAXADD
        DO IA=1,NLIMI
          IF (IA.EQ.NSSPL(N1ST+N2ND+N3RD+IRA)) THEN
            RNUMB(N1ST+N2ND+N3RD+IA)=PRMSPL(N1ST+N2ND+N3RD+IRA)
            NLSPLT(N1ST+N2ND+N3RD+IA)=.TRUE.
          ENDIF
        ENDDO
421   CONTINUE
C
500   IF (.NOT.TRCNAL) RETURN
      WRITE (iunout,*) 
     .  'SPLITTING SURFACES, BELONGING TO THE STANDARD  MESH '
      CALL LEER(1)
      WRITE (iunout,*) 'NLSPLT(ISURF)=TRUE INDICATES THAT THE'
      WRITE (iunout,*) 'SURFACE WITH NUMBER ISURF IS SPLITTING'
      IF (NR1ST.GT.1) THEN
        WRITE (iunout,*) 'RADIAL SURFACES: '
        CALL MASAL1 ('NLSPLT',NLSPLT(1),MIN(NR1ST,N1ST+N2ND+N3RD+NLIM))
      ENDIF
      IF (NP2ND.GT.1) THEN
        WRITE (iunout,*) 'POLOIDAL SURFACES: '
        CALL MASAL1 ('NLSPLT',NLSPLT(N1ST+1),NP2ND)
      ENDIF
      IF (NT3RD.GT.1) THEN
        WRITE (iunout,*) 'TOROIDAL SURFACES: '
        CALL MASAL1 ('NLSPLT',NLSPLT(N1ST+N2ND+1),NT3RD)
      ENDIF
      IF (NLIMI.GE.1) THEN
        WRITE (iunout,*) 'ADDITIONAL SURFACES: '
        CALL MASAL1 ('NLSPLT',NLSPLT(N1ST+N2ND+N3RD+1),NLIMI)
      ENDIF
      CALL LEER(1)
      RETURN
      END
C ===== SOURCE: pedist.f
C
      SUBROUTINE PEDIST (XTIM)
C
C   SUBROUTINE PEDIST CALCULATES THE DISTRIBUTION OF PROCESSORS TO
C   STRATA IN CASE THAT THERE ARE MORE PE'S THAN STRATA
C   DISTRIBUTION OF PE'S IS DONE ACCORDING TO THE DISTRIBUTION OF
C   COMPUTATION TIME.
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CPES
      USE COMSOU
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XTIM(0:NSTRA)
      REAL(DP) :: TIMPE(0:NSTRA)
      REAL(DP) :: FACP, DELT, SUMTIM, TMEAN
      INTEGER :: IPE, K, I, ISTRA, NPRS_FREE, NPRS_OPT

! calculate mean cpu time per stratum
      sumtim=xtim(nstrai)
      TMEAN=SUMTIM/FLOAT(NPRS)

      WRITE (iunout,*) ' SUMTIM = ',SUMTIM,' MEAN TIME = ',TMEAN

      NPRS_OPT=0
      NPRS_FREE=NPRS
      DO ISTRA=1,NSTRAI
        delt=xtim(istra)-xtim(istra-1)
        IF (delt.GE.1.E-5) THEN
! a stratum that has got computation time gets at least 1 processor
          NPESTR(ISTRA)=1
          NPRS_FREE=NPRS_FREE-1
        ELSE
          NPESTR(ISTRA)=0
        ENDIF
! calculate the optimal number of additional processors according to
! distribution of cpu time done in mcarlo (according to number of particles
! and source strength specified in the input)
        TIMPE(ISTRA)=MAX(delt-TMEAN,1.E-5_DP)/TMEAN
        NPRS_OPT=NPRS_OPT+int(TIMPE(ISTRA))
      ENDDO
      WRITE (iunout,*) ' ISTRA, TIMPE '
      DO ISTRA=1,NSTRAI
        WRITE (iunout,*) ISTRA,TIMPE(ISTRA)
      ENDDO

      WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE

! distribute free processors to strata by their optimal number of processors
      FACP=MIN(1.D0,REAL(NPRS_FREE,KIND(1.D0))/
     .             (REAL(NPRS_OPT,KIND(1.D0))+eps30))
      write (iunout,*) ' facp ',facp
      NPESTR(0)=NPRS
      DO ISTRA=1,NSTRAI
        NPESTR(ISTRA)=NPESTR(ISTRA)+int(TIMPE(ISTRA)*FACP)
        NPRS_FREE=NPRS_FREE-int(TIMPE(ISTRA)*FACP)
      ENDDO
      WRITE (iunout,*) ' NPESTR ',(NPESTR(ISTRA),ISTRA=1,NSTRAI)
      WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE

! if there are still free processors left distribute them to all
! strata with more than tmean cpu time assigned to them using a
! daisy chain mechanism
      ISTRA=0
      DO WHILE (NPRS_FREE.GT.0)
        ISTRA=ISTRA+1
        IF (ISTRA.GT.NSTRAI) ISTRA=1
        IF (TIMPE(ISTRA).GT.1.E-10) THEN
          NPESTR(ISTRA)=NPESTR(ISTRA)+1
          NPRS_FREE=NPRS_FREE-1
        ENDIF
      ENDDO
      WRITE (iunout,*) ' NPESTR '
      WRITE (iunout,'(12I6)') (NPESTR(ISTRA),ISTRA=1,NSTRAI)
      WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE

! assign each processor the number of the stratum it shall work on
      IPE=0
      DO ISTRA=1,NSTRAI
        DO K=1,NPESTR(ISTRA)
          IPE=IPE+1
          NSTRPE(IPE-1)=ISTRA
        ENDDO
      ENDDO
      WRITE (iunout,*) ' IPE, ISTRA '
      WRITE (iunout,'(12I6)') (I,NSTRPE(I),I=0,NPRS-1)

! for each stratum define the number of the first processor
! that does calculations for this stratum
! this is used to determine the groups of processors in the
! accumulation of the results for one stratum
      NPESTA(0)=0
      NPESTA(1)=0
      DO ISTRA=2,NSTRAI
        NPESTA(ISTRA)=NPESTA(ISTRA-1)+NPESTR(ISTRA-1)
      ENDDO
      WRITE (iunout,*) ' NPESTA '
      WRITE (iunout,'(12I6)') (NPESTA(I),I=0,NSTRAI)

      RETURN
      END
C ===== SOURCE: plasma_deriv.f
c  new in 2004:
c  density models to contruct backgound data from other given data :
c      Saha, Boltzmann, Corona, Colrad, File (fort.13,or: fort.10)
c
c  presently:  "File" and "Boltzmann": may affect electron density.
c              hence: done prior to electron density, etc...
c              "Corona", "Colrad", "Saha": need electron density as
c                            input, or, at least, do not affect n_e
c                            hence: done after electron density, etc...
C  may05
c  1) additional density model only for neutrals?  removed
c  2) boltzmann factor only if Ti gt tvac
c     to be checked: correct low T limit: everything in lower level?
c  3) new ti only if nlmlti
c  4) new vi only if nlmlv (nlmlv is new, in cinit.f, and set in input.f)
c  june05  
c     bug fix: in colrad model: vnew=max(vvac,vold) nonsense, because
c          vold<0 possible. replaced by vnew=vold
c 
c  march 06
c     new option: icall > 0, and call base_density
c       allows to use output tallies and special "density model" to
c       construct new input tallies (densities, temperatures, drift velocities)
c       e.g. for post processing (diagno), or for iterations.
c
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced

c
      SUBROUTINE PLASMA_DERIV (ICALL)
c  input:
c    icall=0
c      called prior to Monte Carlo Loop
c      in this call all density models referring to output tallies
c      are ignored.
c    icall=1
c      called after to Monte Carlo Loop and sum over strata
c        this allows to put output tallies from a run onto the
c        background for a next iteration or post processing.
c        in this call all density models referring to input tallies  are
c        ignored.
c      write fort.13 after all density models are done.
c  set derived plasma parameters:
c   DEIN             : electron density (from quasineutrality)
c   DEINL            : log electron density (with cutoffs)
c   TEINL            : log electron temperature (with cutoffs)
c   LGVAC(...,NPLS+1): electron vacuum flag
c   LGVAC(...,IPLS)  : bulk ion vacuum flag
c   LGVAC(...,0)     : background vacuum flag
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
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
      USE CESTIM
      use csdvi
      use csdvi_bgk
      use csdvi_cop
      use comsou
      use cspei

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICALL
      REAL(DP) :: ZTII, ZTNI, FCT2, FCRG, FCT1, VDION, ZTEI, ZTNE,
     .            EMPLS, FCT0, TEPLS, DEPLS, DIPLS, AM1, TEF, DEF,
     .            TEI, DEJ, TVACL, DVACL, BOLTZFAC, RCORONA, RCOLRAD,
     .            TEIDEJ, RATE_COEFF, RCMIN, RCMAX
      REAL(DP) :: tpb1, tpb2, second_own
      REAL(DP) :: COEF(0:8), COEF2D(0:8,0:8), FP(6)
      REAL(DP), ALLOCATABLE :: DEINTF(:), SUMNI(:), SUMMNI(:),
     .                         BASE_DENSITY(:), BASE_TEMP(:)
      INTEGER :: I, IR, IN, IP, IPM, J, IPLS, IOLD, ISW, IRE, I1,
     .           IO, IPLSTI, IPLSV, IOLDTI, IOLDV, IBS, JFEXMN, JFEXMX

      TYPE(EIRENE_SPECTRUM), POINTER :: SPEC
      LOGICAL :: FOUND

      FP = 0._DP
      RCMIN = -HUGE(1._DP)
      RCMAX =  HUGE(1._DP)
      JFEXMN = 0
      JFEXMX = 0
        
      tpb1 = second_own()

      IBS = 0
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
          IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMUSR, IO= ',IO
          CLOSE (UNIT=13)
          IF (IO.EQ.0) THEN
            IOLD=TDMPAR(IPLS)%TDM%ISP(1)
            IOLDTI=MPLSTI(IOLD)
            IOLDV=MPLSV(IOLD)
            TIIN(IPLSTI,:)=TIINTF(IOLDTI,:)
            DIIN(IPLS,:)=DIINTF(IOLD,:)
            VXIN(IPLSV,:)=VXINTF(IOLDV,:)
            VYIN(IPLSV,:)=VYINTF(IOLDV,:)
            VZIN(IPLSV,:)=VZINTF(IOLDV,:)
          ENDIF
          DEALLOCATE(DEINTF)
        ELSEIF (INDEX(CDENMODEL(IPLS),'FORT.10') > 0) THEN
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)

          ALLOCATE (BASE_DENSITY(NRAD))
          ALLOCATE (BASE_TEMP(NRAD))
          CALL GET_BASE_DENSITY(1)

          TIIN(IPLSTI,:)=BASE_TEMP(:)
          VXIN(IPLSV,:)=VXIN(IOLDV,:)
          VYIN(IPLSV,:)=VYIN(IOLDV,:)
          VZIN(IPLSV,:)=VZIN(IOLDV,:)
 
          DO IR=1,NSBOX 
            DIIN(IPLS,IR)=MAX(DVAC,BASE_DENSITY(IR))
          ENDDO
          DEALLOCATE (BASE_DENSITY)
          DEALLOCATE (BASE_TEMP)

          IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
            FOUND = .TRUE.
            DO IR=1,NSBOX 
              IF (LSPCCLL(IR)) THEN
                IF (FOUND) ALLOCATE (SPEC)
                CALL GET_SPECTRUM (IR,1,SPEC,FOUND)
                IF (FOUND) THEN
                  IBS = IBS + 1
                  SPEC%IPRTYP = 4
                  SPEC%IPRSP = IPLS
                  BACK_SPEC(IBS)%PSPC => SPEC
                END IF 
              END IF    
            ENDDO
            IF (.NOT.FOUND) DEALLOCATE(SPEC)
          END IF

        ELSEIF (INDEX(CDENMODEL(IPLS),'BOLTZMANN') > 0) THEN
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)

          ALLOCATE (BASE_DENSITY(NRAD))
          ALLOCATE (BASE_TEMP(NRAD))
          CALL GET_BASE_DENSITY(1)

          TIIN(IPLSTI,:)=BASE_TEMP(:)
          VXIN(IPLSV,:)=VXIN(IOLDV,:)
          VYIN(IPLSV,:)=VYIN(IOLDV,:)
          VZIN(IPLSV,:)=VZIN(IOLDV,:)
 
          FOUND = .TRUE.
          DO IR=1,NSBOX 
            BOLTZFAC=0._DP 
            IF (TIIN(IOLDTI,IR).GT.TVAC)
     .        BOLTZFAC=TDMPAR(IPLS)%TDM%G_BOLTZ*
     .                 EXP(-TDMPAR(IPLS)%TDM%DELTAE/BASE_TEMP(IR))
            DIIN(IPLS,IR)=MAX(DVAC,BASE_DENSITY(IR))*BOLTZFAC

            IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
              IF (LSPCCLL(IR)) THEN
                IF (FOUND) ALLOCATE (SPEC)
                CALL GET_SPECTRUM (IR,1,SPEC,FOUND)
                IF (FOUND) THEN
                  IBS = IBS + 1
                  SPEC%SPC = SPEC%SPC * BOLTZFAC
                  SPEC%IPRTYP = 4
                  SPEC%IPRSP = IPLS
                  BACK_SPEC(IBS)%PSPC => SPEC
                END IF 
              END IF    
            END IF

          ENDDO
          IF (.NOT.FOUND) DEALLOCATE(SPEC)
          DEALLOCATE (BASE_DENSITY)
          DEALLOCATE (BASE_TEMP)
        END IF
      END DO

cdr   tpb2 = second_own()
cdr   write (6,*) ' cputime for fort.13, fort.10, boltzmann ',tpb2-tpb1
cdr   tpb1 = tpb2

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
        TEPLS=TEIN(J)
        DEPLS=DEIN(J)
        LGVAC(J,NPLS+1)=TEPLS.LE.TVAC.OR.DEPLS.LE.DVAC
        LGVAC(J,0)     =LGVAC(J,0).AND.LGVAC(J,NPLS+1)
5102  CONTINUE

cdr      tpb2 = second_own()
cdr      write (6,*) ' cputime for log values ',tpb2-tpb1
cdr      tpb1 = tpb2

C
C   SPECIAL DENSITY MODELS
C
      ALLOCATE (BASE_DENSITY(NRAD))
      ALLOCATE (BASE_TEMP(NRAD))
      DO IPLS=1,NPLSI

        IF (LEN_TRIM(CDENMODEL(IPLS)) == 0) CYCLE

        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)

        SELECT CASE (CDENMODEL(IPLS))

        CASE ('SAHA      ')
!PB   TO BE WRITTEN
          WRITE (iunout,*) ' DENSITY PROFILE ACCORDING TO SAHA IS NOT ',
     .                'AVAILABLE '
          WRITE (iunout,*) ' PLEASE CHOOSE DIFFERENT OPTION FOR ION ',
     .                'DENSITY ',IPLS
          CALL EXIT_OWN(1)

        CASE ('CORONA    ')
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)

          CALL GET_BASE_DENSITY(1)

          TIIN(IPLSTI,:)=BASE_TEMP(:)
          VXIN(IPLSV,:)=VXIN(IOLDV,:)
          VYIN(IPLSV,:)=VYIN(IOLDV,:)
          VZIN(IPLSV,:)=VZIN(IOLDV,:)

          REACDAT(NREACI+1)%LRTC = .FALSE.
          CALL SLREAC (NREACI+1,TDMPAR(IPLS)%TDM%FNAME(1),
     .                 TDMPAR(IPLS)%TDM%H2(1),
     .                 TDMPAR(IPLS)%TDM%REACTION(1),
     .                 TDMPAR(IPLS)%TDM%CR(1),
     .                 RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
          AM1=1._DP/TDMPAR(IPLS)%TDM%A_CORONA
c  now COEF contains the excitation rate, 
c  and AMI is the radiative decay rate
c  compute new density of species ipls from equilibrium between
c  these two processes  for the "ground state" density BASE_DENSITY
          FOUND = .TRUE.
          DO IR=1,NSURF
            IF (LGVAC(IR,NPLS+1)) THEN
              TEF=TVACL
            ELSE
              TEF=TEINL(IR)
            END IF
            RCORONA = RATE_COEFF(NREACI+1,TEF,0._DP,.TRUE.,0)
            DIIN(IPLS,IR)=BASE_DENSITY(IR)*RCORONA*DEIN(IR)*AM1
            DIIN(IPLS,IR)=MAX(DVAC,DIIN(IPLS,IR))

            IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
              IF (LSPCCLL(IR)) THEN
                IF (FOUND) ALLOCATE (SPEC)
                CALL GET_SPECTRUM (IR,1,SPEC,FOUND)
                IF (FOUND) THEN
                  IBS = IBS + 1
                  SPEC%IPRTYP = 4
                  SPEC%IPRSP = IPLS
                  SPEC%SPC = SPEC%SPC * RCORONA*DEIN(IR)*AM1
                  BACK_SPEC(IBS)%PSPC => SPEC
                END IF 
              END IF    
            END IF

          END DO
          IF (.NOT.FOUND) DEALLOCATE(SPEC)

cdr      tpb2 = second_own()
cdr      write (6,*) ' cputime for corona ',ipls,tpb2-tpb1
cdr      tpb1 = tpb2

        CASE ('COLRAD    ')
          IF (.NOT.ALLOCATED(SUMNI)) THEN
            ALLOCATE (SUMNI(NRAD))
            ALLOCATE (SUMMNI(NRAD))
          END IF
          SUMNI = 0._DP
          SUMMNI = 0._DP
          DIIN(IPLS,:)=0._DP
C  ARE THERE MULTIPLE ION TEMPERATURES?
          if (nlmlti) then
            TIIN(IPLSTI,:)=0._DP
          endif
C  ARE THERE MULTIPLE ION DRIFT VELOCITIES?
          if (nlmlv) then
            VXIN(IPLSV,:)=0._DP
            VYIN(IPLSV,:)=0._DP
            VZIN(IPLSV,:)=0._DP
          endif
          FOUND = .TRUE.
          DO IRE=1,TDMPAR(IPLS)%TDM%NRE
            IOLD=TDMPAR(IPLS)%TDM%ISP(IRE)
            IOLDTI=MPLSTI(IOLD)
            IOLDV=MPLSV(IOLD)
            CALL GET_BASE_DENSITY(IRE)
            REACDAT(NREACI+1)%LOTH = .FALSE.
            CALL SLREAC (NREACI+1,TDMPAR(IPLS)%TDM%FNAME(IRE),
     .                   TDMPAR(IPLS)%TDM%H2(IRE),
     .                   TDMPAR(IPLS)%TDM%REACTION(IRE),
     .                   TDMPAR(IPLS)%TDM%CR(IRE),
     .                   RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
            I1=INDEX(TDMPAR(IPLS)%TDM%H2(IRE),'.')
            READ (TDMPAR(IPLS)%TDM%H2(IRE)(I1+1:),*) ISW
            SELECT CASE (ISW)
            CASE (11)
c  only temperature dependence in reduced population coefficient
              COEF(0:8)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(1:9,1)
              DO IR=1,NSURF
                IF (LGVAC(IR,NPLS+1)) CYCLE
                TEF=TEINL(IR)
                RCOLRAD=0._DP
                DO I=0,8
                  TEI=TEF**I
                  RCOLRAD=RCOLRAD+COEF(I)*TEI
                END DO
                RCOLRAD=EXP(RCOLRAD)
                DIIN(IPLS,IR)=DIIN(IPLS,IR)+BASE_DENSITY(IR)*RCOLRAD
                if (nlmlti) then
                  TIIN(IPLSTI,IR)=TIIN(IPLSTI,IR)+
     .                            BASE_DENSITY(IR)*BASE_TEMP(IR)
                endif
                if (nlmlv) then
                  VXIN(IPLSV,IR)=VXIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VXIN(IOLDV,IR)
                  VYIN(IPLSV,IR)=VYIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VYIN(IOLDV,IR)
                  VZIN(IPLSV,IR)=VZIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VZIN(IOLDV,IR)
                endif
                SUMNI(IR) = SUMNI(IR) + BASE_DENSITY(IR)
                SUMMNI(IR) = SUMMNI(IR) + NMASSP(IOLD)*BASE_DENSITY(IR)

                IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
                  IF (LSPCCLL(IR)) THEN
                    IF (FOUND) ALLOCATE (SPEC)
                    CALL GET_SPECTRUM (IR,IRE,SPEC,FOUND)
                    IF (FOUND) THEN
                      IF (IRE == 1) THEN
                        IBS = IBS + 1
                        SPEC%IPRTYP = 4
                        SPEC%IPRSP = IPLS
                        SPEC%SPC = SPEC%SPC * RCOLRAD
                        BACK_SPEC(IBS)%PSPC => SPEC
                      ELSE
                        IF (  (BACK_SPEC(IBS)%PSPC%NSPC == SPEC%NSPC)
     .                   .AND.(BACK_SPEC(IBS)%PSPC%SPCMIN==SPEC%SPCMIN) 
     .                   .AND.(BACK_SPEC(IBS)%PSPC%SPCMAX==SPEC%SPCMAX))
     .                  THEN 
                          SPEC%SPC = SPEC%SPC * RCOLRAD
                          BACK_SPEC(IBS)%PSPC%SPC = 
     .                         BACK_SPEC(IBS)%PSPC%SPC + SPEC%SPC
                          DEALLOCATE(SPEC)
                        ELSE
                          WRITE (IUNOUT,*) ' ERROR IN PLASMA_DERIV,',
     .                       ' DENSITY MODEL COLRAD '
                          WRITE (IUNOUT,*) ' TWO SPECTRA CONTRIBTING ',
     .                       ' TO THE SAME BACKGROUND SPECTRUM DO NOT',
     .                       ' MATCH '
                          WRITE (IUNOUT,*) ' IPLS, IRE ',IPLS, IRE 
                          CALL EXIT_OWN(1)
                        END IF
                      END IF
                    END IF 
                  END IF    
                END IF

              END DO
              IF (.NOT.FOUND) DEALLOCATE(SPEC)
            CASE (12)
c  temperature and density dependence in reduced population coefficient
              COEF2D(0:8,0:8)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(1:9,1:9)
              DO IR=1,NSURF
                IF (LGVAC(IR,NPLS+1)) CYCLE
                TEF=TEINL(IR)
                DEF=LOG(DEIN(IR)*1.E-8_DP)
                RCOLRAD=0._DP
                DEJ=1._DP
                DO J=0,8
                  TEIDEJ=DEJ
                  DO I=0,8
                    RCOLRAD=RCOLRAD+COEF2D(I,J)*TEIDEJ
                    TEIDEJ=TEIDEJ*TEF
                  END DO
                  DEJ=DEJ*DEF
                END DO
                RCOLRAD=EXP(RCOLRAD)
                DIIN(IPLS,IR)=DIIN(IPLS,IR)+BASE_DENSITY(IR)*RCOLRAD
                if (nlmlti) then
                  TIIN(IPLSTI,IR)=TIIN(IPLSTI,IR)+
     .                            BASE_DENSITY(IR)*BASE_TEMP(IR)
                endif
                if (nlmlv) then
                  VXIN(IPLSV,IR)=VXIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VXIN(IOLDV,IR)
                  VYIN(IPLSV,IR)=VYIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VYIN(IOLDV,IR)
                  VZIN(IPLSV,IR)=VZIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VZIN(IOLDV,IR)
                endif
                SUMNI(IR) = SUMNI(IR) + BASE_DENSITY(IR)
                SUMMNI(IR) = SUMMNI(IR) + NMASSP(IOLD)*BASE_DENSITY(IR)

                IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
                  IF (LSPCCLL(IR)) THEN
                    IF (FOUND) ALLOCATE (SPEC)
                    CALL GET_SPECTRUM (IR,IRE,SPEC,FOUND)
                    IF (FOUND) THEN
                      IF (IRE == 1) THEN
                        IBS = IBS + 1
                        SPEC%IPRTYP = 4
                        SPEC%IPRSP = IPLS
                        SPEC%SPC = SPEC%SPC * RCOLRAD
                        BACK_SPEC(IBS)%PSPC => SPEC
                      ELSE
                        IF (  (BACK_SPEC(IBS)%PSPC%NSPC == SPEC%NSPC)
     .                   .AND.(BACK_SPEC(IBS)%PSPC%SPCMIN==SPEC%SPCMIN) 
     .                   .AND.(BACK_SPEC(IBS)%PSPC%SPCMAX==SPEC%SPCMAX))
     .                  THEN 
                          SPEC%SPC = SPEC%SPC * RCOLRAD
                          BACK_SPEC(IBS)%PSPC%SPC = 
     .                         BACK_SPEC(IBS)%PSPC%SPC + SPEC%SPC
                          DEALLOCATE(SPEC)
                        ELSE
                          WRITE (IUNOUT,*) ' ERROR IN PLASMA_DERIV,',
     .                       ' DENSITY MODEL COLRAD '
                          WRITE (IUNOUT,*) ' TWO SPECTRA CONTRIBTING ',
     .                       ' TO THE SAME BACKGROUND SPECTRUM DO NOT',
     .                       ' MATCH '
                          WRITE (IUNOUT,*) ' IPLS, IRE ',IPLS, IRE 
                          CALL EXIT_OWN(1)
                        END IF
                      END IF
                    END IF 
                  END IF    
                END IF

              END DO  ! IR
              IF (.NOT.FOUND) DEALLOCATE(SPEC)
            CASE DEFAULT
              WRITE (iunout,*) ' H.',ISW,
     .             ' NOT FORESEEN IN COLRAD DENSITY MODEL '
            END SELECT
          END DO   ! IRE

          DIIN(IPLS,:)=MAX(DVAC,DIIN(IPLS,:))
          if (nlmlti) then
            TIIN(IPLSTI,:)=MAX(TVAC,TIIN(IPLSTI,:)/(SUMNI(:)+eps60))
          endif
          if (nlmlv) then
            VXIN(IPLSV,:)=VXIN(IPLSV,:)/(SUMMNI(:)+eps60)
            VYIN(IPLSV,:)=VYIN(IPLSV,:)/(SUMMNI(:)+eps60)
            VZIN(IPLSV,:)=VZIN(IPLSV,:)/(SUMMNI(:)+eps60)
          endif

        CASE DEFAULT
!  NOTHING TO BE DONE HERE, ALREADY COMPLETED
        END SELECT ! density model

cdr      tpb2 = second_own()
cdr      write (6,*) ' cputime for colrad ',ipls,tpb2-tpb1
cdr      tpb1 = tpb2

      END DO
      IF (ALLOCATED(SUMNI)) THEN
        DEALLOCATE (SUMNI)
        DEALLOCATE (SUMMNI)
      END IF
      DEALLOCATE (BASE_DENSITY)
      DEALLOCATE (BASE_TEMP)

      NBACK_SPEC = IBS

cdr      tpb2 = second_own()
cdr      write (6,*) ' cputime for density models ',tpb2-tpb1
cdr      tpb1 = tpb2

C
C  special density models done

C  SET DRIFT ENERGY (EV)
      DO J=1,NSBOX
        DO IPLS=1,NPLSI
          IPLSV=MPLSV(IPLS)
          IF (NLDRFT) THEN
            IF (INDPRO(4) == 8) THEN
              IF(IPLS.EQ.1) THEN
                EDRIFT(IPLS,J)=CVRSSP(IPLS)*VDION(J)**2
              ELSE
                WRITE(iunout,*)'WARNING! IPLS>1 NO DRIFT!'
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
          DIPLS=DIIN(IPLS,J)
          LGVAC(J,IPLS)=EMPLS.LE.TVAC.OR.DIPLS.LE.DVAC
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
        BVIN(IPLSV,:)=0._DP
        PARMOM(IPLS,:)=0._DP
        DO 5105 J=1,NSBOX
          ZTII=MAX(TVAC,MIN(TIIN(IPLSTI,J),1.E10_DP))
          TIINL(IPLSTI,J)=LOG(ZTII)
          IF (NLDRFT) THEN
            BVIN(IPLSV,J)=BXIN(J)*VXIN(IPLSV,J)+
     .                    BYIN(J)*VYIN(IPLSV,J)+
     .                    BZIN(J)*VZIN(IPLSV,J)
            PARMOM(IPLS,J)=BVIN(IPLSV,J)*SIGN(1._DP,BVIN(IPLSV,J))*
     .                     AMUA*RMASSP(IPLS)
          ENDIF
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

!  

C
C
C  SAVE PLASMA DATA AND ATOMIC DATA ON FORT.13
C

cdr      tpb2 = second_own()
cdr      write (6,*) ' cputime for edrift, b_perp, etc. ',tpb2-tpb1
cdr      tpb1 = tpb2

      IF ((NFILEL >=1) .AND. (NFILEL <=5)) THEN
         NFILEL=3
         CALL WRPLAM(TRCFLE,0)
!      ELSE
      ELSE IF (NFILEL > 5) THEN
         NFILEL=9
         CALL WRPLAM_XDR(TRCFLE,0)
      END IF

cdr      tpb2 = second_own()
cdr      write (6,*) ' cputime for wrplam ',tpb2-tpb1
cdr      tpb1 = tpb2


      RETURN

      CONTAINS

      SUBROUTINE GET_BASE_DENSITY(IRE)
c  input: ire, number of reaction/density that contributes to the
c              evaluation of the expression for the selected species
c              ipls with special density/temperature option
      INTEGER, INTENT(IN) :: IRE
      INTEGER :: IG, IT, ISTRA, ITYP

      BASE_DENSITY = 0._DP
      BASE_TEMP = 0._DP

      IF (ICALL > 0) THEN

C FOR CALLS AFTER PARTICLE TRACING

        ISTRA = TDMPAR(IPLS)%TDM%ISTR(IRE)
        ITYP  = TDMPAR(IPLS)%TDM%ITP(IRE)
        IF (ISTRA.EQ.IESTR.OR.ITYP.EQ.4) THEN
C  NOTHING TO BE DONE
        ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
          IESTR=ISTRA
          CALL RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .               ESTIMV,ESTIMS,ESTIML,
     .               NSDVI1,SDVI1,NSDVI2,SDVI2,
     .               NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .               NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .               NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .               NSIGI_SPC,TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL SYMET(ESTIMV,NTALV,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .                 NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
          IESTR=ISTRA
          CALL RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .               ESTIMV,ESTIMS,ESTIML,
     .               NSDVI1,SDVI1,NSDVI2,SDVI2,
     .               NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .               NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .               NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .               NSIGI_SPC,TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL SYMET(ESTIMV,NTALV,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .                 NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ELSE
          WRITE (6,*) 
     .      'ERROR IN PLASMA_DERIV: DATA FOR STRATUM ISTRA= ',ISTRA
          WRITE (6,*) 'ARE NOT AVAILABLE. '
          RETURN    
        ENDIF
      END IF   ! ICALL > 0

      SELECT CASE (TDMPAR(IPLS)%TDM%ITP(IRE))
      CASE(0)
        IF (ASSOCIATED(PDENPH)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0)  THEN
              BASE_DENSITY(IG) = PDENPH(IOLD,IT)
              BASE_TEMP(IG) = EDENPH(IOLD,IT)/(PDENPH(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF 
      CASE(1)
        IF (ASSOCIATED(PDENA)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0) THEN
              BASE_DENSITY(IG) = PDENA(IOLD,IT)
              BASE_TEMP(IG) = EDENA(IOLD,IT)/(PDENA(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(2)
        IF (ASSOCIATED(PDENM)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0)  THEN
              BASE_DENSITY(IG) = PDENM(IOLD,IT)
              BASE_TEMP(IG) = EDENM(IOLD,IT)/(PDENM(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(3)
        IF (ASSOCIATED(PDENI)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0)  THEN
              BASE_DENSITY(IG) = PDENI(IOLD,IT)
              BASE_TEMP(IG) = EDENI(IOLD,IT)/(PDENI(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF 
      CASE(4)
        BASE_DENSITY(1:NRAD) = DIIN(IOLD,1:NRAD)
        BASE_TEMP(1:NRAD) = TIIN(IOLDTI,1:NRAD)
      CASE DEFAULT
        WRITE (IUNOUT,*) 'BULK SPECIES ', IPLS,
     .                   ' DENSITY MODEL: ',CDENMODEL(IPLS)
        WRITE (IUNOUT,*) ' WRONG PARTICLE TYPE SPECIFIED '
        WRITE (IUNOUT,*) ' CHECK INPUT AND RERUN CASE '
        CALL EXIT_OWN(1)
      END SELECT
         
      END SUBROUTINE GET_BASE_DENSITY



      SUBROUTINE GET_SPECTRUM (ICELL,IRE,SPEC,FOUND)

      INTEGER, INTENT(IN) :: ICELL, IRE
      TYPE(EIRENE_SPECTRUM), INTENT(OUT) :: SPEC
      LOGICAL, INTENT(OUT) :: FOUND
      INTEGER :: ISPC

      FOUND = .FALSE.
      
      DO ISPC = 1, NADSPC
        IF ((ESTIML(ISPC)%PSPC%ISRFCLL == 2) .AND. 
     .      (ESTIML(ISPC)%PSPC%ISPCSRF == ICELL) .AND.
     .      (ESTIML(ISPC)%PSPC%IPRTYP == TDMPAR(IPLS)%TDM%ITP(IRE)).AND.
     .      (ESTIML(ISPC)%PSPC%IPRSP == TDMPAR(IPLS)%TDM%ISP(IRE))) THEN
          ALLOCATE (SPEC%SPC(0:ESTIML(ISPC)%PSPC%NSPC+1))
          ALLOCATE (SPEC%SDV(0:ESTIML(ISPC)%PSPC%NSPC+1))
          ALLOCATE (SPEC%SGM(0:ESTIML(ISPC)%PSPC%NSPC+1))
          SPEC = ESTIML(ISPC)%PSPC
          FOUND = .TRUE.
        END IF
      END DO
       
      END SUBROUTINE GET_SPECTRUM


      END SUBROUTINE PLASMA_DERIV

C ===== SOURCE: plasma.f
C
C
!  6.12.05  bugfix: avoid calculation of B-field in dead cells
!                   because geometrical parameters may not be known there.
!  6.8. 06  bugfix of bugfix: avoid calculation of B-field in dead cells, but
!                             still make sure to set B-field in 1D cases.
      SUBROUTINE PLASMA

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTEXT
      USE COMPRT

      IMPLICIT NONE

      REAL(DP), ALLOCATABLE :: HELP(:)
      REAL(DP) :: PUX, PUY, EL, EP, PN, BD, B, BVAC, FACT
      INTEGER :: IB, IAIN, K, JJ, ITALI, ICELL, IND, ISTREAM,
     .           IP, IT, IA, J, IR, IPLSTI, IPLSV
C
C  INDPRO=9 MEANS: THESE ARRAYS ARE ALREADY SET IN COUPLE_... (SUBR. INFCOP)
      IF (INDPRO(1) /= 9) TEIN = 0.D0
      IF (INDPRO(2) /= 9) TIIN = 0.D0
      DEIN = 0.D0
      IF (INDPRO(3) /= 9) DIIN = 0.D0
      IF (INDPRO(4) /= 9) VXIN = 0.D0
      IF (INDPRO(4) /= 9) VYIN = 0.D0
      IF (INDPRO(4) /= 9) VZIN = 0.D0
      IF (INDPRO(5) /= 9) BXIN = 0.D0
      IF (INDPRO(5) /= 9) BYIN = 0.D0
      IF (INDPRO(5) /= 9) BZIN = 0.D0
      IF (INDPRO(5) /= 9) BFIN = 0.D0
      IF (INDPRO(6) /= 9) ADIN = 0.D0

      ALLOCATE (HELP(NRAD))
      HELP=0.
C
C  SET EIRENE VACUUM MODEL DATA. I.E. IF ALL TEMPERATURES ARE
C  LESS THAN TVAC OR THE BACKGROUND DENSITY IS LESS THAN DVAC,
C  THEN THIS ZONE IS CONSIDERED TO BE AN "EIRENE VACUUM ZONE",
C  PARTICLE MEAN FREE PATHES IN SUCH ZONES ARE
C  EQUAL TO 1.D10 (CM) AND REACTION RATES ARE ZERO
      TVAC=0.02
      DVAC=1.D2
      VVAC=0.
      BVAC=1.
C
C     SET DENSITY, TEMPERATURE AND MACH NUMBER PROFILES
C     ON MESH "RHOZNE(J)"
C
C
C
C  ELECTRON TEMPERATURE
      IND=INDPRO(1)
      GOTO (101,102,103,104,105,106,107,110,110),IND
101     CALL PROFN (TEIN,TE0,TE1,TE2,TE3,TE4,TE5,TVAC)
        GOTO 110
102     CALL PROFE (TEIN,TE0,TE1,TE2,TE4,TE5,TVAC)
        GOTO 110
103     CALL PROFS (TEIN,TE0,TE1,TE5,TVAC)
        GOTO 110
104     CONTINUE
        ISTREAM=TE0
        ITALI=1
        CALL READTL(TXTPLS(1,ITALI),TXTPSP(1,ITALI),
     .              TXTPUN(1,ITALI),
     .              TEIN,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .              3,ISTREAM)
        GOTO 110
105     CALL PROUSR (TEIN,0,TE0,TE1,TE2,TE3,TE4,TE5,TVAC,NSURF)
        GOTO 110
106     CALL PROFR (TEIN,0,1,1,NSURF)
        GOTO 110
107     CALL PROFR (TEIN,0,1,1,NSBOX)
        GOTO 110
110   CONTINUE



C  ION TEMPERATURE
      IND=INDPRO(2)
      DO 120 IPLS=1,NPLSTI
        GOTO (111,112,113,114,115,116,117,120,120),IND
111     CALL PROFN (HELP,TI0(IPLS),TI1(IPLS),TI2(IPLS),
     .                   TI3(IPLS),TI4(IPLS),TI5(IPLS),TVAC)
        TIIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 120
112     CALL PROFE (HELP,TI0(IPLS),TI1(IPLS),TI2(IPLS),
     .                   TI4(IPLS),TI5(IPLS),TVAC)
        TIIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 120
113     CALL PROFS (HELP,TI0(IPLS),TI1(IPLS),TI5(IPLS),TVAC)
        TIIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 120
c  INDPRO=4:  read from stream TIO(IPLS)
114     CONTINUE
        ISTREAM=TI0(IPLS)
        ITALI=2
        CALL READTL(TXTPLS(IPLS,ITALI),TXTPSP(IPLS,ITALI),
     .              TXTPUN(IPLS,ITALI),
     .              HELP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .              3,ISTREAM)
        TIIN(IPLS,1:NSBOX)=HELP(1:NSBOX)
        GOTO 120
115     CALL PROUSR (HELP,1+0*NPLS,TI0(IPLS),TI1(IPLS),TI2(IPLS),
     .                      TI3(IPLS),TI4(IPLS),TI5(IPLS),TVAC,NSURF)
        TIIN(IPLS,1:NSURF)=HELP(1:NSURF)
        GOTO 120
120     CONTINUE
        GOTO 1120
116     CALL PROFR (TIIN,1+0*NPLS,NPLSTI,NPLSTI,NSURF)
        GOTO 1120
117     CALL PROFR (TIIN,1+0*NPLS,NPLSTI,NPLSTI,NSBOX)
        GOTO 1120
1120  CONTINUE



C  ION DENSITY
      IND=INDPRO(3)
      DO 130 IPLS=1,NPLSI
        IF (LEN_TRIM(CDENMODEL(IPLS)) > 0) CYCLE
        GOTO (121,122,123,124,125,126,127,130,130),IND
121     CALL PROFN (HELP,DI0(IPLS),DI1(IPLS),DI2(IPLS),
     .                   DI3(IPLS),DI4(IPLS),DI5(IPLS),DVAC)
        DIIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 130
122     CALL PROFE (HELP,DI0(IPLS),DI1(IPLS),DI2(IPLS),
     .                   DI4(IPLS),DI5(IPLS),DVAC)
        DIIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 130
123     CALL PROFS (HELP,DI0(IPLS),DI1(IPLS),DI5(IPLS),DVAC)
        DIIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 130
c  INDPRO=4:  read from stream DIO(IPLS)
124     CONTINUE
        ISTREAM=DI0(IPLS)
        ITALI=4
        CALL READTL(TXTPLS(IPLS,ITALI),TXTPSP(IPLS,ITALI),
     .              TXTPUN(IPLS,ITALI),
     .              HELP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .              3,ISTREAM)
        DIIN(IPLS,1:NSBOX)=HELP(1:NSBOX)
        GOTO 130
125     CALL PROUSR (HELP,1+1*NPLS,DI0(IPLS),DI1(IPLS),DI2(IPLS),
     .                    DI3(IPLS),DI4(IPLS),DI5(IPLS),DVAC,NSURF)
        DIIN(IPLS,1:NSURF)=HELP(1:NSURF)
        GOTO 130
130     CONTINUE
        GOTO 1130
126     CALL PROFR (DIIN,1+0*NPLS+NPLSTI,NPLSI,NPLS,NSURF)
        GOTO 1130
127     CALL PROFR (DIIN,1+0*NPLS+NPLSTI,NPLSI,NPLS,NSBOX)
        GOTO 1130
1130  CONTINUE



C  DRIFT VELOCITY
      IND=INDPRO(4)
      DO 140 IPLS=1,NPLSV
        GOTO (131,132,133,134,135,136,137,140,140),IND
131     CALL PROFN (HELP,VX0(IPLS),VX1(IPLS),VX2(IPLS),
     .                   VX3(IPLS),VX4(IPLS),VX5(IPLS),VVAC)
        VXIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        CALL PROFN (HELP,VY0(IPLS),VY1(IPLS),VY2(IPLS),
     .                   VY3(IPLS),VY4(IPLS),VY5(IPLS),VVAC)
        VYIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        CALL PROFN (HELP,VZ0(IPLS),VZ1(IPLS),VZ2(IPLS),
     .                   VZ3(IPLS),VZ4(IPLS),VZ5(IPLS),VVAC)
        VZIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 140
132     CALL PROFE (HELP,VX0(IPLS),VX1(IPLS),VX2(IPLS),
     .                   VX4(IPLS),VX5(IPLS),VVAC)
        VXIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        CALL PROFE (HELP,VY0(IPLS),VY1(IPLS),VY2(IPLS),
     .                   VY4(IPLS),VY5(IPLS),VVAC)
        VYIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        CALL PROFE (HELP,VZ0(IPLS),VZ1(IPLS),VZ2(IPLS),
     .                   VZ4(IPLS),VZ5(IPLS),VVAC)
        VZIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 140
133     CALL PROFS (HELP,VX0(IPLS),VX1(IPLS),VX5(IPLS),VVAC)
        VXIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        CALL PROFS (HELP,VY0(IPLS),VY1(IPLS),VY5(IPLS),VVAC)
        VYIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        CALL PROFS (HELP,VZ0(IPLS),VZ1(IPLS),VZ5(IPLS),VVAC)
        VZIN(IPLS,1:NR1ST)=HELP(1:NR1ST)
        GOTO 140
C  INDPRO=4 NOT IN USE
134     CONTINUE
        GOTO 140
135     CALL PROUSR (HELP,1+2*NPLS,VX0(IPLS),VX1(IPLS),VX2(IPLS),
     .                        VX3(IPLS),VX4(IPLS),VX5(IPLS),VVAC,NSURF)
        VXIN(IPLS,1:NSURF)=HELP(1:NSURF)
        CALL PROUSR (HELP,1+3*NPLS,VY0(IPLS),VY1(IPLS),VY2(IPLS),
     .                        VY3(IPLS),VY4(IPLS),VY5(IPLS),VVAC,NSURF)
        VYIN(IPLS,1:NSURF)=HELP(1:NSURF)
        CALL PROUSR (HELP,1+4*NPLS,VZ0(IPLS),VZ1(IPLS),VZ2(IPLS),
     .                        VZ3(IPLS),VZ4(IPLS),VZ5(IPLS),VVAC,NSURF)
        VZIN(IPLS,1:NSURF)=HELP(1:NSURF)
        GOTO 140
140   CONTINUE
C  SCALE FROM MACH NUMBER PROFILE TO CM/SEC PROFILE?
      IF (NLMACH) THEN
        DO 1141 IPLS=1,NPLSI
          IPLSTI=MPLSTI(IPLS)
          IPLSV=MPLSV(IPLS)
          DO 1142 ICELL=1,NSURF
            FACT=CVEL2A*SQRT((TIIN(IPLSTI,ICELL)+
     .                        TEIN(ICELL))/RMASSP(IPLS))
            VXIN(IPLSV,ICELL)=VXIN(IPLSV,ICELL)*FACT
            VYIN(IPLSV,ICELL)=VYIN(IPLSV,ICELL)*FACT
            VZIN(IPLSV,ICELL)=VZIN(IPLSV,ICELL)*FACT
1142      CONTINUE
1141    CONTINUE
      ENDIF
      GOTO 1140
136   CALL PROFR (VXIN,1+1*NPLS+NPLSTI+0*NPLSV,NPLSV,NPLSV,NSURF)
      CALL PROFR (VYIN,1+1*NPLS+NPLSTI+1*NPLSV,NPLSV,NPLSV,NSURF)
      CALL PROFR (VZIN,1+1*NPLS+NPLSTI+2*NPLSV,NPLSV,NPLSV,NSURF)
      GOTO 1140
137   CALL PROFR (VXIN,1+1*NPLS+NPLSTI+0*NPLSV,NPLSV,NPLSV,NSBOX)
      CALL PROFR (VYIN,1+1*NPLS+NPLSTI+1*NPLSV,NPLSV,NPLSV,NSBOX)
      CALL PROFR (VZIN,1+1*NPLS+NPLSTI+2*NPLSV,NPLSV,NPLSV,NSBOX)
      GOTO 1140
1140  CONTINUE
C
C
C  MAGNETIC FIELD UNIT VECTOR
C  FOR IND=5,6,7 OR 9: ALSO THE ABSOLUTE B-FIELD STRENGTH BF CAN BE SET
      IND=INDPRO(5)
C  DEFAULT: BFIELD IN Z-DIRECTION, IE., PITCH=0
      DO J=1,NSBOX
        BXIN(J)=0.
        BYIN(J)=0.
        BZIN(J)=1.
      END DO
      GOTO (141,142,143,144,145,146,147,150,150),IND
C  HELP IS FIELD LINE PITCH ANGLE: B_POL/B_TOT
141     CALL PROFN (HELP,B0,B1,B2,B3,B4,B5,BVAC)
        GOTO 1400
142     CALL PROFE (HELP,B0,B1,B2,B4,B5,BVAC)
        GOTO 1400
143     CALL PROFS (HELP,B0,B1,B5,BVAC)
        GOTO 1400
C  INDPRO=4: read from stream B0:  NOT IN USE
144     CONTINUE
        GOTO 150
C  CONVERT PITCH ANGLE INTO B-FIELD UNIT VECTOR
1400    CONTINUE
        IF (LEVGEO.EQ.1) THEN
          DO 1401 J=1,NSURF
            CALL NCELLN(J,IR,IP,IT,IA,IB,
     .                  NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
            IF (IR.GE.NR1ST) GOTO 1401
            IF ((NP2ND.GT.1).AND.(IP.GE.NP2ND)) GOTO 1401
            BXIN(J)=0.0
            BYIN(J)=HELP(IR)
            BZIN(J)=SQRT(1.-HELP(IR)*HELP(IR))
1401      CONTINUE
        ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
          DO 1402 J=1,NSURF
            CALL NCELLN(J,IR,IP,IT,IA,IB,
     .                  NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
            IF (IR.GE.NR1ST) GOTO 1402
            IF ((NP2ND.GT.1).AND.(IP.GE.NP2ND)) GOTO 1402
            EP=0.5*(EP1(IR)+EP1(IR+1))
            EL=0.5*(ELL(IR)+ELL(IR+1))
            JJ=IR+(IP-1)*NR1ST
            PUY= XCOM(JJ)-EP
            PUX=-YCOM(JJ)/EL/EL
            PN=SQRT(PUX*PUX+PUY*PUY+EPS60)
            BXIN(J)=HELP(IR)*PUX/PN
            BYIN(J)=HELP(IR)*PUY/PN
            BZIN(J)=SQRT(1.-HELP(IR)*HELP(IR))
1402      CONTINUE
        ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
          DO 1403 J=1,NSURF
            CALL NCELLN(J,IR,IP,IT,IA,IB,
     .                  NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
            IF (IR.GE.NR1ST) GOTO 1403
            IF ((NP2ND.GT.1).AND.(IP.GE.NP2ND)) GOTO 1403
            PUX=VPLX(IR,IP)+VPLX(IR+1,IP)
            PUY=VPLY(IR,IP)+VPLY(IR+1,IP)
            PN=SQRT(PUX*PUX+PUY*PUY+EPS60)
            BXIN(J)=HELP(IR)*PUX/PN
            BYIN(J)=HELP(IR)*PUY/PN
            BZIN(J)=SQRT(1.-HELP(IR)*HELP(IR))
1403      CONTINUE
        ELSE
          CALL LEER(1)
          WRITE (iunout,*) 
     .      'DEFAULT MAGNETIC FIELD (IN Z-DIRECTION) IS USED '
          CALL LEER(1)
        ENDIF
        GOTO 150
c  INDPRO=5:  call prousr
145     CONTINUE
        CALL PROUSR (BXIN,1+5*NPLS,B0,B1,B2,B3,B4,B5,0._DP,NSURF)
        CALL PROUSR (BYIN,2+5*NPLS,B0,B1,B2,B3,B4,B5,0._DP,NSURF)
        CALL PROUSR (BZIN,3+5*NPLS,B0,B1,B2,B3,B4,B5,1._DP,NSURF)
        CALL PROUSR (BFIN,4+5*NPLS,B0,B1,B2,B3,B4,B5,1._DP,NSURF)
        GOTO 150
c  INDPRO=6:  call profr (information comes from interfacing code)
c             default (vaccum) parameters in additional cells
146     CALL PROFR (BXIN,1+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        CALL PROFR (BYIN,2+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        CALL PROFR (BZIN,3+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        CALL PROFR (BFIN,4+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
        GOTO 150
c  INDPRO=7:  call profr (information comes from interfacing code)
c             include also additional cells
147     CALL PROFR (BXIN,1+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        CALL PROFR (BYIN,2+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        CALL PROFR (BZIN,3+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        CALL PROFR (BFIN,4+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        GOTO 150
150   CONTINUE
C
C  CHECK FOR ZERO MAGNETIC FIELD IN ANY CELL (INCL. ADD. CELL REGION)
      DO 153 JJ=1,NSBOX
        IF (BXIN(JJ)**2+BYIN(JJ)**2+BZIN(JJ)**2.LE.EPS30) THEN
          WRITE (iunout,*) 'ZERO B-FIELD IN STANDARD CELL JJ= ',JJ
          CALL EXIT_OWN(1)
        ENDIF
        B=SQRT(BXIN(JJ)**2+BYIN(JJ)**2+BZIN(JJ)**2)
        IF (ABS(B-1.D0).GT.EPS12) THEN
          WRITE (iunout,*) 'B-FIELD IN STANDARD CELL JJ= ',JJ,B
          CALL EXIT_OWN(1)
        ENDIF
153   CONTINUE
C
C  ADDITIONAL INPUT TALLIES
      IND=INDPRO(6)
      DO 160 K=1,NAINI
        GOTO (151,151,151,151,155,156,157,160,160),IND
C  DEFAULT: ZERO
151     CONTINUE
        DO 1151 J=1,NR1ST
          ADIN(K,J)=0.
1151    CONTINUE
        GOTO 160
155     CALL PROUSR (HELP,6+5*NPLS,BD,BD,BD,BD,BD,BD,0._DP,NSURF)
!pb        CALL RESETP (ADIN,HELP,K,1,NAIN,NSURF)
        ADIN(K,1:NSURF)=HELP(1:NSURF)
        GOTO 160
160   CONTINUE
      GOTO 1160
156   CALL PROFR (ADIN,6+1*NPLS+NPLSTI+3*NPLSV,NAINI,NAIN,NSURF)
      GOTO 1160
157   CALL PROFR (ADIN,6+1*NPLS+NPLSTI+3*NPLSV,NAINI,NAIN,NSBOX)
      GOTO 1160
1160  CONTINUE
C
C
C   SET VACUUM DATA IN ADDITIONAL REGIONS OUTSIDE THE
C   THE STANDARD MESH
C
      IF (INDPRO(1).LE.6 .OR. INDPRO(1).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          TEIN(J)=TVAC
        ENDDO
      ENDIF
      IF (INDPRO(2).LE.6 .OR. INDPRO(2).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 17 IPLS=1,NPLSI
            IPLSTI=MPLSTI(IPLS)
            TIIN(IPLSTI,J)=TVAC
17        CONTINUE
        ENDDO
      ENDIF
      IF (INDPRO(3).LE.6 .OR. INDPRO(3).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 18 IPLS=1,NPLSI
            DIIN(IPLS,J)=DVAC
18        CONTINUE
        ENDDO
      ENDIF
      IF (INDPRO(4).LE.6 .OR. INDPRO(4).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 19 IPLS=1,NPLSI
            IPLSV=MPLSV(IPLS)
            VXIN(IPLSV,J)=VVAC
            VYIN(IPLSV,J)=VVAC
            VZIN(IPLSV,J)=VVAC
19        CONTINUE
        ENDDO
      ENDIF
      IF (INDPRO(5).LE.6 .OR. INDPRO(5).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          BXIN(J)=0.
          BYIN(J)=0.
          BZIN(J)=1.
        ENDDO
      ENDIF
      IF (INDPRO(6).LE.6 .OR. INDPRO(6).EQ.9) THEN
        DO J=NSURF+1,NSURF+NRADD
          DO 20 IAIN=1,NAINI
            ADIN(IAIN,J)=0.
20        CONTINUE
        ENDDO
      ENDIF

      DEALLOCATE(HELP)
C
      RETURN
      END
C ===== SOURCE: profe.f
C
C
C
C
      SUBROUTINE PROFE(PRO,PRO0,RINP,A1,E,SEP,PROVAC)
C
C  EXPONENTIAL PROFILE, PLUS EXPONENTIAL DECAY BEYOND RHOSEP
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGRID
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PRO0, RINP, A1, E, SEP, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      REAL(DP) :: EP1R, ELLR, YR, AREAA, ARCA, RHOINP, ZR,
     .          AR, DIST, PROS, RHO
      INTEGER :: J, JM1, LEARCA, NLOCAL

      RHO=SEP
      RHOINP=RINP
      IF (LEVGEO.LE.2) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
C FIND RADIAL SURFACE LABELING CO-ORDINATE AT "SEP"
C
        IF (LEVGEO.EQ.2) THEN
          JM1=LEARCA(SEP,RSURF,1,NLOCAL,1,'PROFE       ')
          AR=AREAA (SEP,JM1,ARCA,YR,EP1R,ELLR)
          RHO=SQRT(AR*PIAI)
          JM1=LEARCA(RINP,RSURF,1,NLOCAL,1,'PROFE       ')
          AR=AREAA (RINP,JM1,ARCA,YR,EP1R,ELLR)
          RHOINP=SQRT(AR*PIAI)
        ENDIF
      ELSE
        WRITE (iunout,*) 'WARNING: PROFE CALLED WITH LEVGEO.GT.2 '
        WRITE (iunout,*) 'NO PLASMA PARAMETERS RETURNED'
        RETURN
      ENDIF
C
      DO 20 J=1,NLOCAL-1
        IF (RHOZNE(J).LE.RHO) THEN
          ZR=RHOZNE(J)-RHOINP
          PRO(J)=PRO0*EXP(-ZR/A1)
        ELSE
          DIST=RHO-RHOINP
          PROS=PRO0*EXP(-DIST/A1)
          PRO(J)=PROS*EXP((RHO-RHOZNE(J))/E)
        ENDIF
20    CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END
C ===== SOURCE: profn.f
C sept.05:  activate this radial profile type also for levgeo=3
C           input radius SEP is relative to flux suface labelling grid RHOSRF
C           rra gt. raa indicates: one more vaccum cell in radial direction.

C
      SUBROUTINE PROFN(PRO,PRO0,PROS,P,Q,E,SEP,PROVAC)
C
C  PARABOLIC PROFILE, PLUS EXPONENTIAL DECAY BEYOND RHOSEP
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGRID
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PRO0, PROS, P, Q, E, SEP, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      REAL(DP) :: EP1R, ELLR, YR, AREAA, ARCA, RHOSEP, RORP, FACT, ZR,
     .            RDI, AR
      INTEGER  :: J, JM1, LEARCA, NLOCAL

      IF (LEVGEO.LE.3) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
C FIND RADIAL SURFACE LABELING CO-ORDINATE AT "SEP"
C
        IF (LEVGEO.EQ.2) THEN
          JM1=LEARCA(SEP,RSURF,1,NLOCAL,1,'PROFN       ')
          AR=AREAA (SEP,JM1,ARCA,YR,EP1R,ELLR)
          RHOSEP=SQRT(AR*PIAI)
        ELSEIF (LEVGEO.EQ.1) THEN
          RHOSEP=SEP
        ELSEIF (LEVGEO.EQ.3) THEN
          RHOSEP=SEP
        ENDIF
      ELSE
        WRITE (iunout,*) 'WARNING: SUBR. PROFN CALLED WITH LEVGEO.GT.3 '
        WRITE (iunout,*) 'NO PLASMA PARAMETERS RETURNED'
        RETURN
      ENDIF
C
      RDI=1./(RHOSEP-RHOSRF(1))
      DO 20 J=1,NLOCAL-1
        IF (RHOZNE(J).LE.RHOSEP) THEN
          ZR=RHOZNE(J)-RHOSRF(1)
          RORP=ZR*RDI
          IF (RORP.LE.0.D0) THEN
            PRO(J)=PRO0
          ELSE
            FACT=(1.0-RORP**P)**Q
            PRO(J)=(PRO0-PROS)*FACT+PROS
          ENDIF
        ELSE
          PRO(J)=PROS*EXP((RHOSEP-RHOZNE(J))/E)
        ENDIF
20    CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END
C ===== SOURCE: profr.f
C
C
      SUBROUTINE PROFR (PRO,IINDEX,NSPZI,NSPZ1,NDAT)
C
C  READ ENTIRE PROFILE FROM WORK ARRAY
C
      USE PRECISION
      USE PARMMOD
      USE CSPEI

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IINDEX, NSPZI, NSPZ1, NDAT
      REAL(DP), INTENT(OUT) :: PRO(NSPZ1,*)

      PRO(1:NSPZI,1:NDAT) = PLASMA_BCKGRND(IINDEX+1:IINDEX+NSPZI,1:NDAT)
      RETURN
      END
C ===== SOURCE: profs.f
C
      SUBROUTINE PROFS(PRO,PRO0,PRO1,SEP,PROVAC)
C
C  CONSTANT PROFILE, PLUS SECOND VALUE BEYOND RHOSEP
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGRID
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PRO0, PRO1, SEP, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      REAL(DP) :: AREAA, ARCA, AR, ELLR, YR, EP1R, RHOSEP
      INTEGER :: LEARCA, J, NLOCAL, JM1

      RHOSEP=SEP
      IF (LEVGEO.LE.3) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
C FIND RADIAL SURFACE LABELING CO-ORDINATE AT "SEP"
C
        IF (LEVGEO.EQ.2) THEN
          JM1=LEARCA(SEP,RSURF,1,NLOCAL,1,'PROFS       ')
          AR=AREAA (SEP,JM1,ARCA,YR,EP1R,ELLR)
          RHOSEP=SQRT(AR*PIAI)
        ELSEIF (LEVGEO.EQ.1) THEN
          RHOSEP=SEP
        ELSEIF (LEVGEO.EQ.3) THEN
          RHOSEP=SEP
        ENDIF
      ELSE
        WRITE (iunout,*) 'WARNING: PROFS CALLED WITH LEVGEO.GT.3 '
        WRITE (iunout,*) 'CONSTANT PLASMA PARAMETERS RETURNED'
        NLOCAL=NR1ST
        DO 25 J=1,NLOCAL-1
          PRO(J)=PRO0
25      CONTINUE
        PRO(NR1ST)=0.
        RETURN
      ENDIF
C
      DO 20 J=1,NLOCAL-1
        IF (RHOZNE(J).GT.RHOSEP) GOTO 15
        PRO(J)=PRO0
        GOTO 20
15      PRO(J)=PRO1
20    CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END
C ===== SOURCE: read_tetra.f
      SUBROUTINE READ_TETRA (CASENAME)

!pb 05.12.06: structur coortet is build up from tetrahedra
!pb 07.12.06: set itethand to default value 1

      USE PRECISION
      USE PARMMOD
      USE CTETRA
      USE CLGIN
      USE CGRID
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      CHARACTER*(*), INTENT(IN) :: CASENAME
      CHARACTER(100) :: FILENAME, ZEILE
      INTEGER :: LL, I, IND, IT, IS, IS1, IER, NRK, ISTS,
     .           ISTMIN, ISTMAX, IC, J, JS, JT, IP1, i1, i2, i3, i4
      INTEGER :: ITSIDE(3,4), IP(3), JP(3)
      TYPE(TET_ELEM), POINTER :: CUR
C
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/

      LL=LEN_TRIM(CASENAME)

      FILENAME=CASENAME(1:LL) // '.npco_char'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      READ (ZEILE,*) NRK

      IF (NRK /= NCOORD) THEN
        WRITE (iunout,*) ' NCOORD IS WRONG IN EIRENE INPUT FILE'
        WRITE (iunout,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        CALL EXIT(1)
      END IF

      DO I=1,NCOORD
        READ(30,*) IND, XTETRA(I), YTETRA(I), ZTETRA(I)
      END DO

      CLOSE (UNIT=30)

      FILENAME=CASENAME(1:LL) // '.elemente'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      READ (ZEILE,*) NTET

      IF (NTET > NR1ST) THEN
        WRITE (iunout,*) ' NR1ST IS WRONG IN EIRENE INPUT FILE'
        WRITE (iunout,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        CALL EXIT(1)
      END IF

      DO I=1,NTET
        READ (30,*) IND, NTECK(1,I), NTECK(2,I), NTECK(3,I), NTECK(4,I)
      END DO

      CLOSE (UNIT=30)

      FILENAME=CASENAME(1:LL) // '.neighbors'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      DO I=1,NTET
        READ (30,*) IND, NTBAR(1,I), NTSEITE(1,I), INMTIT(1,I),
     .                   NTBAR(2,I), NTSEITE(2,I), INMTIT(2,I),
     .                   NTBAR(3,I), NTSEITE(3,I), INMTIT(3,I),
     .                   NTBAR(4,I), NTSEITE(4,I), INMTIT(4,I)
        IF (INMTIT(1,I) /= 0) INMTIT(1,I) = INMTIT(1,I) + NLIM
        IF (INMTIT(2,I) /= 0) INMTIT(2,I) = INMTIT(2,I) + NLIM
        IF (INMTIT(3,I) /= 0) INMTIT(3,I) = INMTIT(3,I) + NLIM
        IF (INMTIT(4,I) /= 0) INMTIT(4,I) = INMTIT(4,I) + NLIM
      END DO

      CLOSE (UNIT=30)

      IER = 0
      IF ((MAXVAL(NTECK(1:4,1:NTET)) > NCOORD) .OR.
     .    (MINVAL(NTECK(1:4,1:NTET)) <= 0 )) THEN
        WRITE (iunout,*) ' WRONG COORDINATE NUMBER IS DEFINITION OF',
     .              ' TRIANGLES FOUND '
        IER = 2
      END IF

      IF (.NOT.ALLOCATED(COORTET)) THEN
        ALLOCATE (COORTET(NCOORD))
        DO I=1,NCOORD
          NULLIFY(COORTET(I)%PTET)
        END DO
      END IF

      DO IT=1,NTET
        DO IS=1,4
          IC = NTECK(IS,IT)
          ALLOCATE (CUR)
          CUR%NOTET = IT
          CUR%NEXT_TET => COORTET(IC)%PTET
          COORTET(IC)%PTET => CUR
        ENDDO
      ENDDO

!for testing: setup connection map and write it onto file

!      ntbar = 0
!      ntseite = 0

!      call suche_nachbarn

!      FILENAME=CASENAME(1:LL) // '.neighbors.out'
!      OPEN (UNIT=39,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

!      write (39,'(i10)') ntet 

!      DO I=1,NTET
!        i1 = 0
!        i2 = 0
!        i3 = 0
!        i4 = 0
!        IF (INMTIT(1,I) /= 0) i1 = INMTIT(1,I) - NLIM
!        IF (INMTIT(2,I) /= 0) i2 = INMTIT(2,I) - NLIM
!        IF (INMTIT(3,I) /= 0) i3 = INMTIT(3,I) - NLIM
!        IF (INMTIT(4,I) /= 0) i4 = INMTIT(4,I) - NLIM
!        write (39,'(13i10)') 
!     .        I, NTBAR(1,I), NTSEITE(1,I), i1,
!     .           NTBAR(2,I), NTSEITE(2,I), i2,
!     .           NTBAR(3,I), NTSEITE(3,I), i3,
!     .           NTBAR(4,I), NTSEITE(4,I), i4
!      END DO
!      close (unit=39)

      DO IT=1,NTET
        DO IS=1,4
          IF (NTBAR(IS,IT).EQ.0.AND.INMTIT(IS,IT).EQ.0) THEN
            WRITE (iunout,*) 
            WRITE (iunout,*) ' ERROR IN READ_TETRA '
            WRITE (iunout,*) ' OPEN SIDE OF TETRAHEDRON ',IT,' SIDE ',IS
            WRITE (iunout,*) ' NTECK ',(NTECK(ITSIDE(J,IS),IT),J=1,3)
            IER = 3
          ELSE IF (NTBAR(IS,IT) /= 0) THEN
            JT = NTBAR(IS,IT)
            JS =  NTSEITE(IS,IT)
            IF ((NTBAR(JS,JT) /= IT) .OR. (NTSEITE(JS,JT) /= IS)) THEN
              WRITE (iunout,*) 
              WRITE (iunout,*) ' ERROR IN READ_TETRA '
              WRITE (iunout,*) ' INCONSISTENCY IN CONNECTION MAP'
              WRITE (iunout,*) ' TETRAHEDRON ',IT,' SIDE ',IS,
     .                         ' HAS NEIGHBOR TET ',JT,' SIDE ',JS
              WRITE (iunout,*) ' BUT '
              WRITE (iunout,*) ' TETRAHEDRON ',JT,' SIDE ',JS,
     .                         ' HAS NEIGHBOR TET',NTBAR(JS,JT),
     .                         ' SIDE ',NTSEITE(JS,JT)
              IER = 4
            END IF
            IP(1)=NTECK(ITSIDE(1,IS),IT)
            IP(2)=NTECK(ITSIDE(2,IS),IT)
            IP(3)=NTECK(ITSIDE(3,IS),IT)
            JP(1)=NTECK(ITSIDE(1,JS),JT)
            JP(2)=NTECK(ITSIDE(2,JS),JT)
            JP(3)=NTECK(ITSIDE(3,JS),JT)
            iloop:do i=1,3
              jloop:do j=i,3
                if (ip(i) == jp(j)) then
                  ip1=jp(j)
                  jp(j)=jp(i)
                  jp(i)=ip1
                  cycle iloop
                endif
              end do jloop
              WRITE (iunout,*) 
              write (iunout,*) ' ERROR IN READ_TETRA '
              WRITE (iunout,*) ' INCONSISTENCY IN CONNECTION MAP AND',
     .                         ' ELEMENT DEFINITION '
              WRITE (iunout,*) ' TETRAHEDRON ',IT,' SIDE ',IS,
     .                         ' CONSISTS OF POINTS ', IP         
              WRITE (iunout,*) ' NEIGHBOR TET ',JT,' SIDE ',JS,
     .                         ' CONSISTS OF POINTS ', JP  
              IER = 5
              exit iloop
            end do iloop
          ENDIF
        ENDDO
      ENDDO

      IF (IER /= 0) CALL EXIT(1)

      itethand = 1

      RETURN
      END SUBROUTINE READ_TETRA
C ===== SOURCE: read_triang.f
      SUBROUTINE READ_TRIANG (CASENAME)

      USE PRECISION
      USE PARMMOD
      USE CTRIG
      USE CLGIN
      USE CGRID
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      CHARACTER*(*), INTENT(IN) :: CASENAME
      CHARACTER(100) :: FILENAME, ZEILE
      INTEGER :: LL, I, IND, IT, IS, IS1, IER, NRK, ISTS,
     .           ISTMIN, ISTMAX, IC

      LL=LEN_TRIM(CASENAME)

      FILENAME=CASENAME(1:LL) // '.npco_char'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      READ (ZEILE,*) NRK

      IF (NRK /= NRKNOT) THEN
        WRITE (iunout,*) ' NRKNOT IS WRONG IN EIRENE INPUT FILE'
        WRITE (iunout,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        CALL EXIT_OWN(1)
      END IF

      DO I=1,NRKNOT
        READ(30,*) IND, XTRIAN(I), YTRIAN(I)
      END DO

!pb      XTRIAN(1:NRKNOT) = XTRIAN(1:NRKNOT) * 100._DP
!pb      YTRIAN(1:NRKNOT) = YTRIAN(1:NRKNOT) * 100._DP

      CLOSE (UNIT=30)

      FILENAME=CASENAME(1:LL) // '.elemente'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      READ (ZEILE,*) NTRII

      IF (NTRII > NR1ST) THEN
        WRITE (iunout,*) ' NR1ST IS WRONG IN EIRENE INPUT FILE'
        WRITE (iunout,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        CALL EXIT_OWN(1)
      END IF

      DO I=1,NTRII
        READ (30,*) IND, NECKE(1,I), NECKE(2,I), NECKE(3,I)
      END DO

      CLOSE (UNIT=30)

      FILENAME=CASENAME(1:LL) // '.neighbors'
      OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      END DO

      DO I=1,NTRII
        READ (30,*) IND, NCHBAR(1,I), NSEITE(1,I), INMTI(1,I),
     .                   NCHBAR(2,I), NSEITE(2,I), INMTI(2,I),
     .                   NCHBAR(3,I), NSEITE(3,I), INMTI(3,I),
     .                   IXTRI(I),    IYTRI(I)
        IF (INMTI(1,I) /= 0) INMTI(1,I) = INMTI(1,I) + NLIM
        IF (INMTI(2,I) /= 0) INMTI(2,I) = INMTI(2,I) + NLIM
        IF (INMTI(3,I) /= 0) INMTI(3,I) = INMTI(3,I) + NLIM
      END DO

      CLOSE (UNIT=30)

      IER = 0
      IF ((MAXVAL(NECKE(1:3,1:NTRII)) > NRKNOT) .OR.
     .    (MINVAL(NECKE(1:3,1:NTRII)) <= 0 )) THEN
        WRITE (iunout,*) ' WRONG COORDINATE NUMBER IS DEFINITION OF',
     .              ' TRIANGLES FOUND '
        IER = 2
      END IF

      DO IT=1,NTRII
        DO IS=1,3
          IF (NCHBAR(IS,IT).EQ.0.AND.INMTI(IS,IT).EQ.0) THEN
            WRITE (iunout,*) ' ERROR IN READ_TRIANG '
            WRITE (iunout,*) ' OPEN SIDE OF TRIANGLE ',IT,' SIDE ',IS
            IS1=IS+1
            IF (IS.EQ.3) IS1=1
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS,IT)),
     .                                    YTRIAN(NECKE(IS,IT))
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS1,IT)),
     .                                    YTRIAN(NECKE(IS1,IT))
            IER = 3
          ENDIF
        ENDDO
      ENDDO

      IF (IER /= 0) CALL EXIT_OWN(1)

      ISTMIN=MINVAL(INMTI(1:3,1:NTRII),MASK=(INMTI(1:3,1:NTRII)/=0))
      ISTMAX=MAXVAL(INMTI(1:3,1:NTRII),MASK=(INMTI(1:3,1:NTRII)/=0))

      IC=0
      DO ISTS=ISTMIN,ISTMAX
        DO IS=1,3
          DO IT=1,NTRII
            IF (INMTI(IS,IT) == ISTS) THEN
              IC=IC+1
              INSPAT(IS,IT)=IC
            END IF
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE READ_TRIANG
C ===== SOURCE: setcon.f
C  06.12.05:  added: au_to_cm2 (was previously defined in fpatha)
C  25.08.06:  use double precision constant for definition of PMASSA
c
      SUBROUTINE SETCON
      USE PRECISION
      USE PARMMOD
      USE CCONA

      IMPLICIT NONE

      REAL(DP) :: AMUAKG

      CALL ALLOC_CCONA

C   CONSTANTS
C  ELECTRON CHARGE [EV]
      ELCHA=1.6022D-19
C  ENERGY CONVERSION
      EV_TO_J=ELCHA
      J_TO_EV=1.D0/EV_TO_J
      J_TO_ERG=1.0D7
      ERG_TO_J=1.D0/J_TO_ERG
      EV_TO_ERG=EV_TO_J*J_TO_ERG
      ERG_TO_EV=1.D0/EV_TO_ERG
      EVKEL=8.6173D-5
      EV2HZ  = 2.41798834D14
C  PLANCK CONSTANT, [ERG S]
      HPLANCK=6.6260755D-27
C  PLANCK CONSTANT, [EV S]
      HPLNK=HPLANCK*ERG_TO_EV
C  SPEED OF LIGHT, [CM/S]
      CLIGHT=2.99792458D+10
c h*c [eV*cm] FOR CONVERSION FROM ENERGY TO WAVELENGTH UNITS
      HPCL   = CLIGHT*HPLNK
C  BOHR MAGNETON, [EV/TESLA]
      MUB=5.788381804D-5
C  ATOMIC MASS UNIT [G]
      AMUA=1.6606D-24
C  PROTON, ELECTRON MASS
      PMASSA=1.0073D0
      PMASSE=5.448D-4
C  CONVERT CROSS SECTIONS FROM ATOMIC UNITS TO CM**2
      AU_TO_CM2=5.29177E-9**2
C  NUMERICAL PRECISION PARAMETERS
      EPS60=1.D-60
      EPS30=1.D-30
      EPS12=1.D-12
      EPS10=1.D-10
      EPS6=1.D-6
      EPS5=1.D-5
C
      PIA=4._DP*ATAN(1._DP)
      PI2A=2._DP*PIA
      PIHA=PIA/2._DP
      PIQU=PIA/4._DP
      PISQ=SQRT(PIA)
      PIAI=1._DP/PIA
C
      SQ2=SQRT(2._DP)
      SQ2I=1._DP/SQ2
      DEGRAD=PIA/180._DP
      RADDEG=180._DP/PIA
C  EIRENE UNITS CONVERSIONS
      AMUAKG=AMUA*1.D-3
C  V_THERMAL=CVEL2A*SQRT(T(EV)/RMASS(AMU)),  CVEL2A=0.98227E6
      CVEL2A=SQRT(1.D4*ELCHA/AMUAKG)
C  VEL(CM/S)=CVELAA*SQRT(E0(EV)/RMASS(AMU)), CVELAA=1.38912E6
      CVELAA=SQ2*CVEL2A

      CVELI2=1._DP/CVELAA/CVELAA
      EFACT=CVELI2*PMASSA
      EFCT23=EFACT*2._DP/3._DP
      HPLNK_BAR=HPLNK/PI2A
C
C  IONIZATION POTENTIAL OF NEUTRAL HYDROGEN ATOM
      EIONH=13.6_DP
C  IONIZATION POTENTIAL OF NEUTRAL HYDROGEN MOLECULE
      EIONH2=15.4_DP
C  IONIZATION POTENTIAL OF NEUTRAL HELIUM ATOM
      EIONHE=24.588_DP
C
      RETURN
C
      END
C ===== SOURCE: seteq.f
C
C
! 6.12.05  Bugfix: one index I replaced by J. Bug was relevant for switching of 
!                  IGJUM2 in case of bit arrays
 
      SUBROUTINE SETEQ

      USE PRECISION
      USE PARMMOD
      USE COMPRT, ONLY: IUNOUT
      USE COMUSR
      USE CADGEO
      USE CTRCEI
      USE CLGIN

      IMPLICIT NONE

      REAL(DP) :: SG
      INTEGER :: NLBT, IEQ, J, IDIMP, K, I, IEQ1
      LOGICAL BITGET
C
C  FIRST: ADDITIONAL SURFACES
C
      DO 1 J=1,NLIMI
        IEQ=IABS(ILEQUI(J))
        IF (IEQ.EQ.0) GOTO 1
        IF (IGJUM0(J)+IGJUM0(IEQ) .NE. 0) GOTO 1
        SG=ISIGN(1,ILEQUI(J))
C   COEFFICIENTS OF SURFACE J ARE SET EQUAL TO THOSE OF SURFACE IEQ
        A0LM(J)=SG*A0LM(IEQ)
        A1LM(J)=SG*A1LM(IEQ)
        A2LM(J)=SG*A2LM(IEQ)
        A3LM(J)=SG*A3LM(IEQ)
        A4LM(J)=SG*A4LM(IEQ)
        A5LM(J)=SG*A5LM(IEQ)
        A6LM(J)=SG*A6LM(IEQ)
        A7LM(J)=SG*A7LM(IEQ)
        A8LM(J)=SG*A8LM(IEQ)
        A9LM(J)=SG*A9LM(IEQ)
        IF (JUMLIM(J).EQ.0) THEN
          ALM(J)=SG*ALM(IEQ)
          BLM(J)=SG*BLM(IEQ)
          CLM(J)=SG*CLM(IEQ)
        ELSE
          ALM(J)=ALM(IEQ)
          BLM(J)=BLM(IEQ)
          CLM(J)=CLM(IEQ)
        ENDIF
        IF (JUMLIM(J).NE.0) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM1(J,IEQ)=1
            IGJUM1(IEQ,J)=1
          ELSE
            CALL BITSET (IGJUM1,0,NLIMPS,J,IEQ,1,NBITS)
            CALL BITSET (IGJUM1,0,NLIMPS,IEQ,J,1,NBITS)
          END IF
        ELSE
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM2(J,IEQ)=1
            IGJUM2(IEQ,J)=1
          ELSE
            CALL BITSET (IGJUM2,0,NLIMPS,J,IEQ,1,NBITS)
            CALL BITSET (IGJUM2,0,NLIMPS,IEQ,J,1,NBITS)
          END IF
        ENDIF
C
        IF (TRCSUR) THEN
          WRITE (iunout,*) 'EQUATION OF SURFACE NO. ',J,
     .                     ' IS SET EQUAL TO '
          WRITE (iunout,*) 'EQUATION OF SURFACE NO. ',IEQ
          CALL LEER(1)
        ENDIF
1     CONTINUE
C
C  NEXT: NON DEFAULT STANDARD SURFACES
C  JUMLIM=0 FOR STANDARD SURFACES, ONLY LGJUM2 IS USED IN TIME-ROUTINES
C
      DO 10 J=NLIM+1,NLIM+NSTSI
        IF (ILEQUI(J).EQ.0.OR.IGJUM0(J).NE.0) GOTO 10
        IF (INUMP(J-NLIM,1).NE.0) IDIMP=1
        IF (INUMP(J-NLIM,2).NE.0) IDIMP=2
        IF (INUMP(J-NLIM,3).NE.0) IDIMP=3
C  NEGATIVE SIGN OF ILEQUI IS MEANINGLESS FOR STANDARD SURFACES
C  BECAUSE THEIR ORIENTATION CANNOT BE CHANGED
        IEQ1=IABS(ILEQUI(J))
C  FIND INDEX IEQ OF SURFACE IEQ1, TO BE EQUALLED WITH SURFACE J
C  SEARCH ONLY IN THE SAME (RADIAL, POLOIDAL OR TOROIDAL) GRID
        IEQ=0
        DO 11 I=1,NSTSI
          IF (INUMP(I,IDIMP).EQ.IEQ1) IEQ=NLIM+I
11      CONTINUE
        IF (NLIMPB >= NLIMPS) THEN
          IGJUM2(J,IEQ)=1
          IGJUM2(IEQ,J)=1
        ELSE
          CALL BITSET (IGJUM2,0,NLIMPS,J,IEQ,1,NBITS)
          CALL BITSET (IGJUM2,0,NLIMPS,IEQ,J,1,NBITS)
        END IF
        IF (TRCSUR) THEN
          WRITE (iunout,*) 'EQUATION OF SURFACE NO. ',J,
     .                     ' IS SET EQUAL TO '
          WRITE (iunout,*) 'EQUATION OF SURFACE NO. ',IEQ
          CALL LEER(1)
        ENDIF
10    CONTINUE
C
C  LGJUM1: COMMUTATIVE
C  LGJUM2: ASSOCIATIVE
      IF (NLIMPB >= NLIMPS) THEN
        DO J=1,NLIMI
          IF (SUM(IGJUM1(J,1:NLIMI)) > 1 ) THEN
            DO I=1,NLIMI
              IF (IGJUM1(J,I) .NE. 0) IGJUM1(I,J)=1
            ENDDO
          ENDIF
          IF (JUMLIM(J).EQ.0) THEN
            IF (SUM(IGJUM2(J,1:NLIMI)) > 1 ) THEN
              DO I=1,NLIMI
                DO K=1,NLIMI
                  IF (IGJUM2(I,J)*IGJUM2(J,K) .NE. 0) THEN
                    IGJUM2(I,K)=1
                    IGJUM2(K,I)=1
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ELSE
        NLBT = NLIMI/NBITS+1
        DO J=1,NLIMI
          IF (SUM(IGJUM1(J,1:NLBT)) > 1 ) THEN
            DO I=1,NLIMI
              IF (BITGET(IGJUM1,0,NLIMPS,J,I,NBITS))
     .          CALL BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
            END DO
          END IF
          IF (JUMLIM(J).EQ.0) THEN
            IF (SUM(IGJUM2(J,1:NLBT)) > 1 ) THEN
              DO I=1,NLIMI
                IF (BITGET(IGJUM2,0,NLIMPS,I,J,NBITS)) THEN
                  DO K=1,NLIMI
                    IF (BITGET(IGJUM1,0,NLIMPS,J,K,NBITS)) THEN
                      CALL BITSET (IGJUM2,0,NLIMPS,I,K,1,NBITS)
                      CALL BITSET (IGJUM2,0,NLIMPS,K,I,1,NBITS)
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END IF
        END DO
      END IF
      DO 200 J=NLIM+1,NLIM+NSTSI
        DO 200 I=NLIM+1,NLIM+NSTSI
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM1(J,I).NE.0) IGJUM1(I,J)=1
          DO K=NLIM+1,NLIM+NSTSI
            IF (IGJUM2(I,J)*IGJUM2(J,K) .NE. 0) THEN
              IGJUM2(I,K)=1
              IGJUM2(K,I)=1
            ENDIF
          END DO
        ELSE
          IF (BITGET(IGJUM1,0,NLIMPS,J,I,NBITS))
     .      CALL BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
          DO K=NLIM+1,NLIM+NSTSI
            IF (BITGET(IGJUM2,0,NLIMPS,I,J,NBITS) .AND.
     .          BITGET(IGJUM1,0,NLIMPS,J,K,NBITS)) THEN
              CALL BITSET (IGJUM2,0,NLIMPS,I,K,1,NBITS)
              CALL BITSET (IGJUM2,0,NLIMPS,K,I,1,NBITS)
            END IF
          END DO
        END IF
200   CONTINUE
      RETURN
      END
C ===== SOURCE: setfit.f
      SUBROUTINE SETFIT(TRCSUR)
c new option, 2002: ilfit.lt.0   not fully tested
c                   connecting surface "ie" not necessarily rlb=1, or 1.5
c                   anymore.
      USE PRECISION
      USE PARMMOD
      USE COMPRT, ONLY: IUNOUT
      USE CADGEO
      USE CCONA
      USE CLGIN

      IMPLICIT NONE
C
      LOGICAL, INTENT(IN) :: TRCSUR
      REAL(DP) :: XS(8), YS(8), D, DST(8,2), XA(14), XB(14)
      REAL(DP) :: XT, SRAD, RAD, YMIN, YMAX, B, AH, YT,
     .          DX1, DY2, DZ3, XNORM, AT, A3, A4, A2, XR, XL, A5,
     .          A1, A0, YR, YL, B0, B1, B2, B3, B4, XMAX, XMIN, XP1,
     .          XP2, XP3, XP4, YP1, YP2, YP3, YP4, B5
      INTEGER :: IEQ(2), ILFT, ILFTS
      INTEGER :: IMIN1, IMIN2, K, IS, JUM, IPNT1, IPNT2, I, IE, J
      LOGICAL :: LINFX,LINFY,LINFZ, TWOPOINT
C
      DO 1 I=1,NLIMI
C  IS ILFIT OPTION IN USE?
        IF (ILFIT(I).EQ.0.OR.IGJUM0(I).NE.0) GOTO 1
C  SELECT THE SURFACE NUMBERS OF THE NEIGHBORING SURFACES
        IF (ILFIT(I).LT.0) THEN
          ILFT=-ILFIT(I)
          ILFTS=-1
        ELSE
          ILFT=ILFIT(I)
          ILFTS=1
        ENDIF
        IEQ(1)=ILFT/1000
        IEQ(2)=ILFT-IEQ(1)*1000
C  SURFACE I MUST BE GIVEN BY TWO POINT OPTION WITH ONE
C  IGNORABLE CO-ORDINATE
C  USE THIRD POINT FOR IDENTIFICATION OF 2-POINT INPUT OPTION
C  BECAUSE RLB HAS ALREADY BEEN OVERWRITTEN
        IF (P3(1,I).LT.1.D50.AND.P3(2,I).LT.1.D50.AND.
     .      P3(3,I).LT.1.D50) GOTO 991
C  WHICH CO-ORDINATE OF SURFACE NO. I IS IGNORABLE?
        LINFX=P3(1,I).GT.1.D50
        LINFY=P3(2,I).GT.1.D50
        LINFZ=P3(3,I).GT.1.D50
        IPNT1=0
        IPNT2=0
        DO 2 J=1,2
          IE=IEQ(J)
          IF (IE.EQ.0) GOTO 2
          IF (IGJUM0(IE).NE.0) GOTO 991
C  CONNECT SURFACE I WITH SURFACE NO. IE
C  SURFACE IE MUST BE GIVEN WITH SAME IGNORABLE CO-ORDINATE AS SURFACE I
C
C  IDENTIFY IGNORABLE CO-ORDINATE OF SURFACE IE
          IF (LINFX) THEN
            IF (A1LM(IE).NE.0..OR.A4LM(IE).NE.0..OR.
     .          A7LM(IE).NE.0..OR.A8LM(IE).NE.0.D0) GOTO 5
C  X IS IGNORABLE, IN BOTH SURFACES: I AND IE
            A0=A0LM(IE)
            A1=A2LM(IE)
            A2=A3LM(IE)
            A3=A5LM(IE)
            A4=A6LM(IE)
            A5=A9LM(IE)
            XL=YLIMS1(1,IE)
            XR=YLIMS2(1,IE)
            YL=ZLIMS1(1,IE)
            YR=ZLIMS2(1,IE)

            XP1=P1(2,I)
            XP2=P2(2,I)
            YP1=P1(3,I)
            YP2=P2(3,I)
C
            XP3=P1(2,IE)
            XP4=P2(2,IE)
            YP3=P1(3,IE)
            YP4=P2(3,IE)
C
            GOTO 100
          ENDIF
5         CONTINUE
          LINFX=.FALSE.
          IF (LINFY) THEN
            IF (A2LM(IE).NE.0..OR.A5LM(IE).NE.0..OR.
     .          A7LM(IE).NE.0..OR.A9LM(IE).NE.0.D0) GOTO 6
C  Y IS IGNORABLE, IN BOTH SURFACES: I AND IE
            A0=A0LM(IE)
            A1=A1LM(IE)
            A2=A3LM(IE)
            A3=A4LM(IE)
            A4=A6LM(IE)
            A5=A8LM(IE)
            XL=XLIMS1(1,IE)
            XR=XLIMS2(1,IE)
            YL=ZLIMS1(1,IE)
            YR=ZLIMS2(1,IE)

            XP1=P1(1,I)
            XP2=P2(1,I)
            YP1=P1(3,I)
            YP2=P2(3,I)
C
            XP3=P1(1,IE)
            XP4=P2(1,IE)
            YP3=P1(3,IE)
            YP4=P2(3,IE)
C
            GOTO 100
          ENDIF
6         CONTINUE
          LINFY=.FALSE.
          IF (LINFZ) THEN
            IF (A3LM(IE).NE.0..OR.A6LM(IE).NE.0..OR.
     .          A8LM(IE).NE.0..OR.A9LM(IE).NE.0.D0) GOTO 7
C  Z IS IGNORABLE, IN BOTH SURFACES: I AND IE
            A0=A0LM(IE)
            A1=A1LM(IE)
            A2=A2LM(IE)
            A3=A4LM(IE)
            A4=A5LM(IE)
            A5=A7LM(IE)
            XL=XLIMS1(1,IE)
            XR=XLIMS2(1,IE)
            YL=YLIMS1(1,IE)
            YR=YLIMS2(1,IE)

            XP1=P1(1,I)
            XP2=P2(1,I)
            YP1=P1(2,I)
            YP2=P2(2,I)
C
            XP3=P1(1,IE)
            XP4=P2(1,IE)
            YP3=P1(2,IE)
            YP4=P2(2,IE)
C
            GOTO 100
          ENDIF
7         CONTINUE
          LINFZ=.FALSE.
          GOTO 990
100       CONTINUE
C
C
C  IS SURFACE IE GIVEN BY TWO-POINT OPTION ?
        TWOPOINT=.FALSE.
        IF (P3(1,IE).GE.1.D50.OR.P3(2,IE).GE.1.D50.OR.
     .      P3(3,IE).GE.1.D50) TWOPOINT=.TRUE.

          IF ((RLB(IE).EQ.1..OR.RLB(IE).EQ.1.5).AND..NOT.TWOPOINT) THEN
C  TRY TO CONNECT SURFACE I TO POINTS ON BOUNDARY BOX OF SURFACE IE
          IS=0
C  SCHNITTPUNKTE VON IE MIT X=XL AND X=XR
          IF (ABS(A4).LE.EPS12) THEN
            YT=-(A0+A1*XL+A3*XL*XL)/(A2+A5*XL)
            IF (YT.GE.YL.AND.YT.LE.YR) THEN
              IS=IS+1
              XS(IS)=XL
              YS(IS)=YT
            ENDIF
            YT=-(A0+A1*XR+A3*XR*XR)/(A2+A5*XR)
            IF (YT.GE.YL.AND.YT.LE.YR) THEN
              IS=IS+1
              XS(IS)=XR
              YS(IS)=YT
            ENDIF
          ELSE
            AH=0.5*(A2+A5*XL)/A4
            B=(A0+A1*XL+A3*XL*XL)/A4
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              YT=-AH+SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XL
                YS(IS)=YT
              ENDIF
              YT=-AH-SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XL
                YS(IS)=YT
              ENDIF
            ENDIF
            AH=0.5*(A2+A5*XR)/A4
            B=(A0+A1*XR+A3*XR*XR)/A4
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              YT=-AH+SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XR
                YS(IS)=YT
              ENDIF
              YT=-AH-SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XR
                YS(IS)=YT
              ENDIF
            ENDIF
          ENDIF
C  SCHNITTPUNKTE MIT Y=YL AND Y=YR
          IF (ABS(A3).LE.EPS12) THEN
            XT=-(A0+A2*YL+A4*YL*YL)/(A1+A5*YL)
            IF (XT.GE.XL.AND.XT.LE.XR) THEN
              IS=IS+1
              YS(IS)=YL
              XS(IS)=XT
            ENDIF
            XT=-(A0+A2*YR+A4*YR*YR)/(A1+A5*YR)
            IF (XT.GE.XL.AND.XT.LE.XR) THEN
              IS=IS+1
              YS(IS)=YR
              XS(IS)=XT
            ENDIF
          ELSE
            AH=0.5*(A1+A5*YL)/A3
            B=(A0+A2*YL+A4*YL*YL)/A3
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              XT=-AH+SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YL
                XS(IS)=XT
              ENDIF
              XT=-AH-SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YL
                XS(IS)=XT
              ENDIF
            ENDIF
            AH=0.5*(A1+A5*YR)/A3
            B=(A0+A2*YR+A4*YR*YR)/A3
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              XT=-AH+SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YR
                XS(IS)=XT
              ENDIF
              XT=-AH-SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YR
                XS(IS)=XT
              ENDIF
            ENDIF
          ENDIF

          ELSEIF ((RLB(IE).EQ.1..OR.RLB(IE).EQ.1.5).AND.TWOPOINT) THEN
C  TRY TO CONNECT SURFACE I TO SURFACE IE
            IS=1
            CALL SCHNITP(XP1,YP1,XP2,YP2,XP3,YP3,XP4,YP4,XS(1),YS(1))
          ELSE
C  ALL OTHER RLB(IE) OPTIONS
            GOTO 991
          ENDIF
C
C  SELECT THE DISTANCES TO THE INTERSECTION POINTS
          DO 3 K=1,IS
            DST(K,1)=1.D60
            DST(K,2)=1.D60
            IF (IPNT1.EQ.0)
     .      DST(K,1)=SQRT((XS(K)-XP1)**2+(YS(K)-YP1)**2)
            IF (IPNT2.EQ.0)
     .      DST(K,2)=SQRT((XS(K)-XP2)**2+(YS(K)-YP2)**2)
3         CONTINUE
C
          IMIN1=1
          IMIN2=1
          DO 4 K=2,IS
            IF (DST(K,1).LT.DST(IMIN1,1)) IMIN1=K
            IF (DST(K,2).LT.DST(IMIN2,2)) IMIN2=K
4         CONTINUE
C
          IF (TRCSUR) THEN
            WRITE (iunout,*) 'SETFIT, I,IE ',I,IE
            DO 4711 K=1,IS
              WRITE (iunout,*) 'K,XS,YS,DIST1,DIST2 ',
     .                     K,XS(K),YS(K),DST(K,1),DST(K,2)
4711        CONTINUE
          ENDIF
C
C  SET THE SELECTED POINT
          IF (ILFTS.LT.0) THEN
C  EXCHANGE POINTS TO BE RESET: NOT THE CLOSEST, BUT THE OTHER ONE
            D=DST(IMIN1,1)
            DST(IMIN1,1)=DST(IMIN2,2)
            DST(IMIN2,2)=D
          ENDIF
          IF (DST(IMIN1,1).LT.DST(IMIN2,2)) THEN
            XP1=XS(IMIN1)
            YP1=YS(IMIN1)
            IPNT1=1
            IF (TRCSUR)
     .      WRITE (iunout,*) 'PNT.1 REPLACED BY ',XP1,YP1,
     .                       ' FOR SURF. ',I
          ELSE
            XP2=XS(IMIN2)
            YP2=YS(IMIN2)
            IPNT2=1
            IF (TRCSUR)
     .      WRITE (iunout,*) 'PNT.2 REPLACED BY ',XP2,YP2,
     .                       ' FOR SURF. ',I
          ENDIF
C
C  RESET THE POINTS ON THE ARRAYS
          IF (TRCSUR) WRITE (iunout,*) ' LINFXYZ ',LINFX,LINFY,LINFZ
          IF (LINFX) THEN
            P1(2,I)=XP1
            P1(3,I)=YP1
            P2(2,I)=XP2
            P2(3,I)=YP2
            DZ3=(P1(3,I)-P2(3,I))
            DY2=(P1(2,I)-P2(2,I))
            A0LM(I)=DZ3*P1(2,I)-DY2*P1(3,I)
            A1LM(I)=0.
            A2LM(I)=-DZ3
            A3LM(I)=DY2
            YLIMS1(1,I)=MIN(P1(2,I),P2(2,I))
            YLIMS2(1,I)=MAX(P1(2,I),P2(2,I))
            ZLIMS1(1,I)=MIN(P1(3,I),P2(3,I))
            ZLIMS2(1,I)=MAX(P1(3,I),P2(3,I))
            IF (P1(2,I).EQ.P2(2,I)) THEN
              YLIMS1(1,I)=YLIMS1(1,I)-0.1
              YLIMS2(1,I)=YLIMS2(1,I)+0.1
            ENDIF
            IF (P1(3,I).EQ.P2(3,I)) THEN
              ZLIMS1(1,I)=ZLIMS1(1,I)-0.1
              ZLIMS2(1,I)=ZLIMS2(1,I)+0.1
            ENDIF
          ELSEIF (LINFY) THEN
            P1(1,I)=XP1
            P1(3,I)=YP1
            P2(1,I)=XP2
            P2(3,I)=YP2
            DZ3=(P1(3,I)-P2(3,I))
            DX1=(P1(1,I)-P2(1,I))
            A0LM(I)=DZ3*P1(1,I)-DX1*P1(3,I)
            A1LM(I)=-DZ3
            A2LM(I)=0.
            A3LM(I)=DX1
            XLIMS1(1,I)=MIN(P1(1,I),P2(1,I))
            XLIMS2(1,I)=MAX(P1(1,I),P2(1,I))
            ZLIMS1(1,I)=MIN(P1(3,I),P2(3,I))
            ZLIMS2(1,I)=MAX(P1(3,I),P2(3,I))
            IF (P1(1,I).EQ.P2(1,I)) THEN
              XLIMS1(1,I)=XLIMS1(1,I)-0.1
              XLIMS2(1,I)=XLIMS2(1,I)+0.1
            ENDIF
            IF (P1(3,I).EQ.P2(3,I)) THEN
              ZLIMS1(1,I)=ZLIMS1(1,I)-0.1
              ZLIMS2(1,I)=ZLIMS2(1,I)+0.1
            ENDIF
          ELSEIF (LINFZ) THEN
            P1(1,I)=XP1
            P1(2,I)=YP1
            P2(1,I)=XP2
            P2(2,I)=YP2
            DY2=(P1(2,I)-P2(2,I))
            DX1=(P1(1,I)-P2(1,I))
            A0LM(I)=DY2*P1(1,I)-DX1*P1(2,I)
            A1LM(I)=-DY2
            A2LM(I)=DX1
            A3LM(I)=0.
            XLIMS1(1,I)=MIN(P1(1,I),P2(1,I))
            XLIMS2(1,I)=MAX(P1(1,I),P2(1,I))
            YLIMS1(1,I)=MIN(P1(2,I),P2(2,I))
            YLIMS2(1,I)=MAX(P1(2,I),P2(2,I))
            IF (P1(1,I).EQ.P2(1,I)) THEN
              XLIMS1(1,I)=XLIMS1(1,I)-0.1
              XLIMS2(1,I)=XLIMS2(1,I)+0.1
            ENDIF
            IF (P1(2,I).EQ.P2(2,I)) THEN
              YLIMS1(1,I)=YLIMS1(1,I)-0.1
              YLIMS2(1,I)=YLIMS2(1,I)+0.1
            ENDIF
          ENDIF
C
          AT=MAX(ABS(A1LM(I)),ABS(A2LM(I)),ABS(A3LM(I)))
          IF (ABS(A1LM(I)).EQ.AT) JUMLIM(I)=1
          IF (ABS(A2LM(I)).EQ.AT) JUMLIM(I)=2
          IF (ABS(A3LM(I)).EQ.AT) JUMLIM(I)=3
          XNORM=SQRT(A1LM(I)*A1LM(I)+A2LM(I)*A2LM(I)+A3LM(I)*A3LM(I))
          IF (XNORM.LE.EPS60) GOTO 993
          A0LM(I)=A0LM(I)/XNORM
          A1LM(I)=A1LM(I)/XNORM
          A2LM(I)=A2LM(I)/XNORM
          A3LM(I)=A3LM(I)/XNORM
          JUM=JUMLIM(I)
          GOTO (91,92,93),JUM
91          ALM(I)=-A0LM(I)/A1LM(I)
            BLM(I)=-A2LM(I)/A1LM(I)
            CLM(I)=-A3LM(I)/A1LM(I)
          GOTO 97
92          ALM(I)=-A0LM(I)/A2LM(I)
            BLM(I)=-A1LM(I)/A2LM(I)
            CLM(I)=-A3LM(I)/A2LM(I)
          GOTO 97
93          ALM(I)=-A0LM(I)/A3LM(I)
            BLM(I)=-A1LM(I)/A3LM(I)
            CLM(I)=-A2LM(I)/A3LM(I)
97        CONTINUE
C
          IF (TRCSUR) THEN
            WRITE (iunout,*) ' A0-A3 ',A0LM(I),A1LM(I),A2LM(I),A3LM(I)
            WRITE (iunout,*) ' XLIMS ',XLIMS1(1,I),XLIMS2(1,I)
            WRITE (iunout,*) ' YLIMS ',YLIMS1(1,I),YLIMS2(1,I)
            WRITE (iunout,*) ' ZLIMS ',ZLIMS1(1,I),ZLIMS2(1,I)
          ENDIF
C
C  TWO POINT OPTION FINISHED. P1,P2 REDEFINED
C  ALL OTHER SURFACE COEFFICIENTS ALSO REDEFINED

C  NOW: GENERAL SECOND ORDER EQUATION, RLB=1., FOR SURFACE I
C
2       CONTINUE
1     CONTINUE
C
      RETURN
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. SETFIT '
      WRITE (iunout,*) 
     .  'INCONSISTENCY IN IGNORABLE CO-ORDINATES DETECTED'
      WRITE (iunout,*) 'BETWEEN REQUESTING SURFACE I= ',I,' AND IE= ',IE
      WRITE (iunout,*) 'JUMLIM(I),LINFX,LINFY,LINFZ ',
     .                  JUMLIM(I),LINFX,LINFY,LINFZ
      WRITE (iunout,*) 'EXIT CALLED '
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. SETFIT '
      WRITE (iunout,*) 'FIT OPTION FOR RLB(IE) = ',RLB(IE),
     .                 ' NOT FORESEEN'
      WRITE (iunout,*) 'REQUEST FROM SURFACE NO. ',I
      WRITE (iunout,*) 'IE = ',IE,' EXIT CALLED '
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. SETFIT '
      WRITE (iunout,*) 'THE VALID AREAS DO NOT INTERSECT'
      WRITE (iunout,*) 'I,IE ',I,IE
      CALL EXIT_OWN(1)
993   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. SETFIT '
      WRITE (iunout,*) 'STRAIGHT LINE NO I= ',I,' COLLAPSED TO A POINT'
      WRITE (iunout,*) 'SURFACE NO. I IS REDUNDANT. USE CH0 I/I OPTION '
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: setprm.f - from Petra, 21.07.08 (see file in this directory)

C
C  *************************
C  *                       *
C  * INITIALIZATION PHASE  *
C  *                       *
C  *************************
C
C      SUBROUTINE SETPRM
C
      SUBROUTINE SETPRM
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CESTIM
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CSDVI
      USE CLGIN
      USE COUTAU
      USE COMXS
      USE CSPEI
      USE CTEXT
      USE CTRCEI

      IMPLICIT NONE

      REAL(DP) :: RSAVE
      INTEGER :: ISAVE, NTESTP, J, NTEST, NTESTI, ITAL, NLSTTL, NLSTTW
      LOGICAL :: LEXTALV(NTALV), LEXTALS(NTALS), 
     .           LEXGENA, LEXGENM, LEXGENI, LEXGENPH
C

      IF (NSTORDT.LT.1 .OR. NSTORDT.GT.9) THEN
        WRITE (iunout,*) 
     .    'POSSIBLE STORAGE CONFLICT. NSTORDT OUT OF RANGE '
        CALL EXIT_OWN(1)
      ENDIF
 
 
C  SWITCH OFF SOME VOLUME AVERAGED OUTPUT TALLIES AUTOMATICALLY; 
C  TRY TO KEEP ONLY THOSE TALLIES THAT ARE NEEDED FOR THE TYPE OF 
C  SPECIES PRESENT IN THE PARTICULAR CASE. 
c   e.g.  no photon tallies unless photons are included (NPHOT>0) 
c   e.g.  no test-ion tallies unless test ions are included (NION>0)   
c   e.g.  no generation limit tallies unless there is, indeed a  
c         generation limit activated   

      LEXGENA  = ANY(NGENA(1:NATM) /= 0)
      LEXGENM  = ANY(NGENM(1:NMOL) /= 0)
      LEXGENI  = ANY(NGENI(1:NION) /= 0)
      LEXGENPH = ANY(NGENPH(1:NPHOT) /= 0)

      LEXTALV(1)  =  NATM>0
      LEXTALV(2)  =  NMOL>0
      LEXTALV(3)  =  NION>0
      LEXTALV(4)  =  NPHOT>0
      LEXTALV(5)  =  NATM>0
      LEXTALV(6)  =  NMOL>0
      LEXTALV(7)  =  NION>0
      LEXTALV(8)  =  NPHOT>0 

      LEXTALV(9)  =  NATM>0
      LEXTALV(10) =  NATM>0
      LEXTALV(11) = (NATM>0).and.(NMOL>0)
      LEXTALV(12) = (NATM>0).and.(NION>0)
      LEXTALV(13) = (NATM>0).and.(NPHOT>0)
      LEXTALV(14) = (NATM>0).and.(NPLS>0) 

      LEXTALV(15) =  NMOL>0
      LEXTALV(16) = (NMOL>0).and.(NATM>0)
      LEXTALV(17) =  NMOL>0
      LEXTALV(18) = (NMOL>0).and.(NION>0)
      LEXTALV(19) = (NMOL>0).and.(NPHOT>0)
      LEXTALV(20) = (NMOL>0).and.(NPLS>0) 

      LEXTALV(21) =  NION>0
      LEXTALV(22) = (NION>0).and.(NATM>0)
      LEXTALV(23) = (NION>0).and.(NMOL>0)
      LEXTALV(24) =  NION>0
      LEXTALV(25) = (NION>0).and.(NPHOT>0)
      LEXTALV(26) = (NION>0).and.(NPLS>0) 

      LEXTALV(27) =  NPHOT>0
      LEXTALV(28) = (NPHOT>0).and.(NATM>0)
      LEXTALV(29) = (NPHOT>0).and.(NMOL>0)
      LEXTALV(30) = (NPHOT>0).and.(NION>0)
      LEXTALV(31) =  NPHOT>0
      LEXTALV(32) = (NPHOT>0).and.(NPLS>0) 

      LEXTALV(33) =                NATM>0
      LEXTALV(34) =                NATM>0
      LEXTALV(35) = (NMOL>0) .and.(NATM>0)
      LEXTALV(36) = (NION>0) .and.(NATM>0)
      LEXTALV(37) = (NPHOT>0).and.(NATM>0)
      LEXTALV(38) = (NPLS>0) .and.(NATM>0) 

      LEXTALV(39) =                NMOL>0
      LEXTALV(40) = (NATM>0) .and.(NMOL>0)
      LEXTALV(41) =                NMOL>0
      LEXTALV(42) = (NION>0) .and.(NMOL>0)
      LEXTALV(43) = (NPHOT>0).and.(NMOL>0)
      LEXTALV(44) = (NPLS>0) .and.(NMOL>0) 

      LEXTALV(45) =                NION>0
      LEXTALV(46) = (NATM>0) .and.(NION>0)
      LEXTALV(47) = (NMOL>0) .and.(NION>0)
      LEXTALV(48) =                NION>0
      LEXTALV(49) = (NPHOT>0).and.(NION>0)
      LEXTALV(50) = (NPLS>0) .and.(NION>0) 

      LEXTALV(51) =               NPHOT>0
      LEXTALV(52) = (NATM>0).and.(NPHOT>0)
      LEXTALV(53) = (NMOL>0).and.(NPHOT>0)
      LEXTALV(54) = (NION>0).and.(NPHOT>0)
      LEXTALV(55) =               NPHOT>0
      LEXTALV(56) = (NPLS<0).and.(NPHOT>0) 

      LEXTALV(NTALA) = NADV>0
      LEXTALV(NTALC) = NCLV>0
      LEXTALV(NTALT) = NSNV>0
      LEXTALV(NTALM) = NCPV>0
      LEXTALV(NTALB) = NBGV>0
      LEXTALV(NTALR) = NALV>0
C  GENERATION LIMIT TALLIES 
C  some of these tallies may be 
c  turned off, depending upon whether generation limits 
c  are activated or not.
      LEXTALV(63) = NATM>0  .AND. LEXGENA
      LEXTALV(64) = NMOL>0  .AND. LEXGENM
      LEXTALV(65) = NION>0  .AND. LEXGENI
      LEXTALV(66) = NPHOT>0 .AND. LEXGENPH
      LEXTALV(67) = NATM>0  .AND. LEXGENA
      LEXTALV(68) = NMOL>0  .AND. LEXGENM
      LEXTALV(69) = NION>0  .AND. LEXGENI
      LEXTALV(70) = NPHOT>0 .AND. LEXGENPH
      LEXTALV(71) = NATM>0  .AND. LEXGENA
      LEXTALV(72) = NMOL>0  .AND. LEXGENM
      LEXTALV(73) = NION>0  .AND. LEXGENI
      LEXTALV(74) = NPHOT>0 .AND. LEXGENPH 
C  PRIMARY SOURCE RATES
      LEXTALV(75) = NATM>0
      LEXTALV(76) = NMOL>0
      LEXTALV(77) = NION>0
      LEXTALV(78) = NPHOT>0
      LEXTALV(79) = NPLS>0 

      LEXTALV(80) = NATM>0
      LEXTALV(81) = NMOL>0
      LEXTALV(82) = NION>0
      LEXTALV(83) = NPHOT>0
      LEXTALV(84) = NPLS>0

      LEXTALV(85) = NATM>0
      LEXTALV(86) = NMOL>0
      LEXTALV(87) = NION>0
      LEXTALV(88) = NPHOT>0
      LEXTALV(89) = NATM>0
      LEXTALV(90) = NMOL>0
      LEXTALV(91) = NION>0
      LEXTALV(92) = NPHOT>0
      LEXTALV(93) = NATM>0
      LEXTALV(94) = NMOL>0
      LEXTALV(95) = NION>0
      LEXTALV(96) = NPHOT>0
      LEXTALV(97) = (NPLS>0) .AND. (NATM>0)
      LEXTALV(98) = (NPLS>0) .AND. (NMOL>0)
      LEXTALV(99) = (NPLS>0) .AND. (NION>0)
      LEXTALV(100) = (NPLS>0) .AND. (NPHOT>0)
      

      LIVTALV = LEXTALV .AND. .NOT.LMISTALV
      LMISTALV = LMISTALV .AND. LEXTALV

      LEA = LEAAT .OR. LEAML .OR. LEAIO .OR. LEAPHT .OR. LEAPL
      LEM = LEMAT .OR. LEMML .OR. LEMIO .OR. LEMPHT .OR. LEMPL
      LEIO = LEIAT .OR. LEIML .OR. LEIIO .OR. LEIPHT .OR. LEIPL
      LEPH = LEPHAT .OR. LEPHML .OR. LEPHIO .OR. LEPHPHT .OR. LEPHPL
 
C 
C  LEADING DIMENSIONS OF FIELDS IN COMMON BLOCK CESTIM AND COUTAU 
C                (i.e. of volume- or surface averaged output tallies) 
C 
C  DENSITIES; ENERGY DENSITIES
      NFIRST(1)=NATM
      NFIRST(2)=NMOL
      NFIRST(3)=NION
      NFIRST(4)=NPHOT
      NFIRST(5)=NATM
      NFIRST(6)=NMOL
      NFIRST(7)=NION
      NFIRST(8)=NPHOT 
C  PARTICLE SOURCES
      NFIRST(9)=0
      NFIRST(10)=NATM
      NFIRST(11)=NMOL
      NFIRST(12)=NION
      NFIRST(13)=NPHOT
      NFIRST(14)=NPLS
      NFIRST(15)=0
      NFIRST(16)=NATM
      NFIRST(17)=NMOL
      NFIRST(18)=NION
      NFIRST(19)=NPHOT
      NFIRST(20)=NPLS
      NFIRST(21)=0
      NFIRST(22)=NATM
      NFIRST(23)=NMOL
      NFIRST(24)=NION
      NFIRST(25)=NPHOT
      NFIRST(26)=NPLS
      NFIRST(27)=0
      NFIRST(28)=NATM
      NFIRST(29)=NMOL
      NFIRST(30)=NION
      NFIRST(31)=NPHOT
      NFIRST(32)=NPLS 
C  ENERGY SOURCES
      NFIRST(33)=0
      NFIRST(34)=0
      NFIRST(35)=0
      NFIRST(36)=0
      NFIRST(37)=0
      NFIRST(38)=0
      NFIRST(39)=0
      NFIRST(40)=0
      NFIRST(41)=0
      NFIRST(42)=0
      NFIRST(43)=0
      NFIRST(44)=0
      NFIRST(45)=0
      NFIRST(46)=0
      NFIRST(47)=0
      NFIRST(48)=0
      NFIRST(49)=0
      NFIRST(50)=0
      NFIRST(51)=0
      NFIRST(52)=0
      NFIRST(53)=0
      NFIRST(54)=0
      NFIRST(55)=0
      NFIRST(56)=0
      NFIRST(NTALA)=NADV
      NFIRST(NTALC)=NCLV
      NFIRST(NTALT)=NSNV
      NFIRST(NTALM)=NCPV
      NFIRST(NTALB)=NBGV
      NFIRST(NTALR)=NALV
C  GENERATION LIMIT TALLIES
      NFIRST(63)=NATM
      NFIRST(64)=NMOL
      NFIRST(65)=NION
      NFIRST(66)=NPHOT
      NFIRST(67)=NATM
      NFIRST(68)=NMOL
      NFIRST(69)=NION
      NFIRST(70)=NPHOT
      NFIRST(71)=NATM
      NFIRST(72)=NMOL
      NFIRST(73)=NION
      NFIRST(74)=NPHOT 
C  PRIMARY SOURCE TALLIES
      NFIRST(75)=NATM
      NFIRST(76)=NMOL
      NFIRST(77)=NION
      NFIRST(78)=NPHOT
      NFIRST(79)=NPLS
      NFIRST(80)=0
      NFIRST(81)=0
      NFIRST(82)=0
      NFIRST(83)=0
      NFIRST(84)=0
      NFIRST(85)=NATM
      NFIRST(86)=NMOL
      NFIRST(87)=NION
      NFIRST(88)=NPHOT 
      NFIRST(89)=NATM
      NFIRST(90)=NMOL
      NFIRST(91)=NION
      NFIRST(92)=NPHOT 
      NFIRST(93)=NATM
      NFIRST(94)=NMOL
      NFIRST(95)=NION
      NFIRST(96)=NPHOT 
      NFIRST(97)=NPLS
      NFIRST(98)=NPLS
      NFIRST(99)=NPLS
      NFIRST(100)=NPLS 
C
C  NTALV=100 ?
C
      DO 1 J=1,NTALV
        NFRSTI(J)=NFIRST(J)+1
        IF (LIVTALV(J)) NFIRST(J)=MAX0(1,NFIRST(J))
1     CONTINUE
C
      NADDI(1)=0
      NADDV(1)=0
      NLSTTL = 0
      DO 2 J=2,NTALV
        NADDI(J)=NADDI(J-1)+NFRSTI(J-1)
        IF (LIVTALV(J-1)) THEN
          NADDV(J)=NADDV(J-1)+NFIRST(J-1)
          NLSTTL = J-1
        ELSE
          NADDV(J)=NADDV(J-1)
        END IF
2     CONTINUE

      IF (LIVTALV(NTALV)) NLSTTL = NTALV
C
C  TOTAL NUMBER OF VOLUME AVERAGED TALLIES
!pb      NVOLTL=NADDV(NTALV)+NFIRST(NTALV)
      NVOLTL=NADDV(NTALV)+NFIRST(NLSTTL)

      NTEST=NADDV(NTALV)+NFIRST(NTALV)
      NTESTI=NADDI(NTALV)+NFRSTI(NTALV)
      NTEST=NTEST*NRTAL
      NTESTI=NTESTI*NSTRAP


      LEXTALS(1) =                (NATM>0)
      LEXTALS(2) =                (NATM>0)
      LEXTALS(3) = (NMOL>0)  .and.(NATM>0)
      LEXTALS(4) = (NION>0)  .and.(NATM>0)
      LEXTALS(5) = (NPHOT>0) .and.(NATM>0)
      LEXTALS(6) = (NPLS>0)  .and.(NATM>0)  

      LEXTALS(7) =                (NMOL>0) 
      LEXTALS(8) = (NATM>0)  .and.(NMOL>0)
      LEXTALS(9) =                (NMOL>0)
      LEXTALS(10) =(NION>0)  .and.(NMOL>0)
      LEXTALS(11) =(NPHOT>0) .and.(NMOL>0)
      LEXTALS(12) =(NPLS>0)  .and.(NMOL>0) 

      LEXTALS(13) =               (NION>0)
      LEXTALS(14) =(NATM>0)  .and.(NION>0)
      LEXTALS(15) =(NMOL>0)  .and.(NION>0)
      LEXTALS(16) =               (NION>0)
      LEXTALS(17) =(NPHOT>0) .and.(NION>0)
      LEXTALS(18) =(NPLS>0)  .and.(NION>0) 

      LEXTALS(19) =               (NPHOT>0)
      LEXTALS(20) =(NATM>0)  .and.(NPHOT>0)
      LEXTALS(21) =(NMOL>0)  .and.(NPHOT>0)
      LEXTALS(22) =(NION>0)  .and.(NPHOT>0)
      LEXTALS(23) =               (NPHOT>0)
      LEXTALS(24) =(NPLS>0)  .and.(NPHOT>0) 

      LEXTALS(25) = NPLS>0

      LEXTALS(26) =               (NATM>0)
      LEXTALS(27) =               (NATM>0)
      LEXTALS(28) =(NMOL>0)  .and.(NATM>0) 
      LEXTALS(29) =(NION>0)  .and.(NATM>0)
      LEXTALS(30) =(NPHOT>0) .and.(NATM>0)
      LEXTALS(31) =(NPLS>0)  .and.(NATM>0) 

      LEXTALS(32) =               (NMOL>0)
      LEXTALS(33) =(NATM>0)  .and.(NMOL>0)
      LEXTALS(34) =               (NMOL>0)
      LEXTALS(35) =(NION>0)  .and.(NMOL>0)
      LEXTALS(36) =(NPHOT>0) .and.(NMOL>0)
      LEXTALS(37) =(NPLS>0)  .and.(NMOL>0) 
 
      LEXTALS(38) =               (NION>0)
      LEXTALS(39) =(NATM>0)  .and.(NION>0)
      LEXTALS(40) =(NMOL>0)  .and.(NION>0)
      LEXTALS(41) =               (NION>0)
      LEXTALS(42) =(NPHOT>0) .and.(NION>0)
      LEXTALS(43) =(NPLS>0)  .and.(NION>0)
      
      LEXTALS(44) =               (NPHOT>0)
      LEXTALS(45) =(NATM>0)  .and.(NPHOT>0)
      LEXTALS(46) =(NMOL>0)  .and.(NPHOT>0)
      LEXTALS(47) =(NION>0)  .and.(NPHOT>0)
      LEXTALS(48) =               (NPHOT>0)
      LEXTALS(49) =(NPLS>0)  .and.(NPHOT>0) 

      LEXTALS(50) = NPLS>0
C  SPUTTERED FLUXES; BY INCIDENT SPECIES
      LEXTALS(51) = NATM>0
      LEXTALS(52) = NMOL>0
      LEXTALS(53) = NION>0
      LEXTALS(54) = NPHOT>0
      LEXTALS(55) = NPLS>0 
C  SPUTTERED FLUX; TOTAL
      LEXTALS(56) = .TRUE. 

      LEXTALS(NTLSA) = NADS>0
      LEXTALS(NTLSR) = NALS>0
      LEXTALS(59) = NSPZ>0

      LIVTALS = LEXTALS .AND. .NOT.LMISTALS
      LMISTALS = LMISTALS .AND. LEXTALS
C
      NFRSTW(1)=NATM
      NFRSTW(2)=NATM
      NFRSTW(3)=NATM
      NFRSTW(4)=NATM
      NFRSTW(5)=NATM
      NFRSTW(6)=NATM
      NFRSTW(7)=NMOL
      NFRSTW(8)=NMOL
      NFRSTW(9)=NMOL
      NFRSTW(10)=NMOL
      NFRSTW(11)=NMOL
      NFRSTW(12)=NMOL
      NFRSTW(13)=NION
      NFRSTW(14)=NION
      NFRSTW(15)=NION
      NFRSTW(16)=NION
      NFRSTW(17)=NION
      NFRSTW(18)=NION
      NFRSTW(19)=NPHOT
      NFRSTW(20)=NPHOT
      NFRSTW(21)=NPHOT
      NFRSTW(22)=NPHOT
      NFRSTW(23)=NPHOT
      NFRSTW(24)=NPHOT
      NFRSTW(25)=NPLS

      NFRSTW(26)=NATM
      NFRSTW(27)=NATM
      NFRSTW(28)=NATM
      NFRSTW(29)=NATM
      NFRSTW(30)=NATM
      NFRSTW(31)=NATM
      NFRSTW(32)=NMOL
      NFRSTW(33)=NMOL
      NFRSTW(34)=NMOL
      NFRSTW(35)=NMOL
      NFRSTW(36)=NMOL
      NFRSTW(37)=NMOL
      NFRSTW(38)=NION
      NFRSTW(39)=NION
      NFRSTW(40)=NION
      NFRSTW(41)=NION
      NFRSTW(42)=NION
      NFRSTW(43)=NION
      NFRSTW(44)=NPHOT
      NFRSTW(45)=NPHOT
      NFRSTW(46)=NPHOT
      NFRSTW(47)=NPHOT
      NFRSTW(48)=NPHOT
      NFRSTW(49)=NPHOT
      NFRSTW(50)=NPLS

      NFRSTW(51)=NATM
      NFRSTW(52)=NMOL
      NFRSTW(53)=NION
      NFRSTW(54)=NPHOT
      NFRSTW(55)=NPLS
      NFRSTW(56)=0
      NFRSTW(NTLSA)=NADS
      NFRSTW(NTLSR)=NALS
      NFRSTW(59)=NSPZ
C
C  NTALS=59 ?
C
      DO 4 J=1,NTALS
        NFRTWI(J)=NFRSTW(J)+1
        IF (LIVTALS(J)) NFRSTW(J)=MAX0(1,NFRSTW(J))
4     CONTINUE
C
      NDDWI(1)=0
      NADDW(1)=0
      NLSTTW = 0
      DO 3 J=2,NTALS
        NDDWI(J)=NDDWI(J-1)+NFRTWI(J-1)
        IF (LIVTALS(J-1)) THEN
          NADDW(J)=NADDW(J-1)+NFRSTW(J-1)
          NLSTTW = J-1
        ELSE
          NADDW(J)=NADDW(J-1)
        END IF
3     CONTINUE

      IF (LIVTALS(NTALS)) NLSTTW = NTALS

C  TOTAL NUMBER OF SURFACE AVERAGED TALLIES
!pb      NSRFTL=NADDW(NTALS)+NFRSTW(NTALS)
      NSRFTL=NADDW(NTALS)+NFRSTW(NLSTTW)

      NTEST=NADDW(NTALS)+NFRSTW(NTALS)
      NTESTI=NDDWI(NTALS)+NFRTWI(NTALS)
      NTEST=NTEST*NLMPGS
      NTESTI=NTESTI*NSTRAP

      CALL ALLOC_CESTIM(2)
      CALL ASSOCIATE_CESTIM

C
C  CHECK LENGTH OF ARRAYS, WHICH ARE EQUIVALENZED TO OTHER ARRAYS
C
C
      IF (.FALSE.) THEN
      RSAVE=SGMS(NSD)
      SGMS(NSD)=1.234567
      IF (SDVI1(NSD,NRTAL+1).NE.1.234567) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: NSDVI1?'
        CALL EXIT_OWN(1)
      ENDIF
      SGMS(NSD)=RSAVE
C
      RSAVE=SGMWS(NSDW)
      SGMWS(NSDW)=1.234567
      IF (SDVI2(NSDW,NLIMPS+1).NE.1.234567) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: NSDVI2?'
        CALL EXIT_OWN(1)
      ENDIF
      SGMWS(NSDW)=RSAVE
C
      RSAVE=VGENPH(NPHOT,NRTAL)
      VGENPH(NPHOT,NRTAL)=1.234567
      IF (ESTIMV(NVOLTL,NRTAL).NE.1.234567) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: NESTM1?'
        CALL EXIT_OWN(1)
      ENDIF
      VGENPH(NPHOT,NRTAL)=RSAVE
C
      RSAVE=SPUMP(NSPZ,NLMPGS)
      SPUMP(NSPZ,NLMPGS)=1.234567
      IF (ESTIMS(NSRFTL,NLMPGS).NE.1.234567) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: NESTM2?'
        CALL EXIT_OWN(1)
      ENDIF
      SPUMP(NSPZ,NLMPGS)=RSAVE
C
C  NOW ATOMIC DATA ARRAYS: COMXS
C
      RSAVE=ZMFPI
      ZMFPI=1.234567
      IF (XSTORV(NSTORV).NE.1.234567) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: NSTOR?'
        CALL EXIT_OWN(1)
      ENDIF
      ZMFPI=RSAVE
C
      RSAVE=VOLTOT
      VOLTOT=1.234567
      IF (RCGM1(NCGM1).NE.1.234567) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: NCGM1?'
        CALL EXIT_OWN(1)
      ENDIF
      VOLTOT=RSAVE
C
      RSAVE=ZDF
      ZDF=1.234567
      IF (RCGRID(NCGRD).NE.1.234567) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: NCGRD?'
        CALL EXIT_OWN(1)
      ENDIF
      ZDF=RSAVE
C
      ISAVE=NSBOX_TAL
      NSBOX_TAL=1234567
      IF (ICGRID(MCGRD).NE.1234567) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: MCGRD?'
        CALL EXIT_OWN(1)
      ENDIF
      NSBOX_TAL=ISAVE
C
      ISAVE=NPPLG
      NPPLG=1234567
      IF (ICPLYG(MCPLYG).NE.1234567) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: MCPLYG?'
        CALL EXIT_OWN(1)
      ENDIF
      NPPLG=ISAVE
      END IF
C
      NFRSTP(1)=0
      NFRSTP(2)=NPLSTI
      NFRSTP(3)=0
      NFRSTP(4)=NPLS
      NFRSTP(5)=NPLSV
      NFRSTP(6)=NPLSV
      NFRSTP(7)=NPLSV
      NFRSTP(8)=0
      NFRSTP(9)=0
      NFRSTP(10)=0
      NFRSTP(11)=0
      NFRSTP(12)=NAIN
      NFRSTP(13)=NPLS
      NFRSTP(14)=0
      NFRSTP(15)=NSPZMC
      NFRSTP(16)=0
      NFRSTP(17)=0
C
C  NTALI=17?
C
      DO 5 J=1,NTALI
        NFRSTP(J)=MAX0(1,NFRSTP(J))
5     CONTINUE
C
      NADDP(1)=0
      DO 6 J=2,NTALI
6       NADDP(J)=NADDP(J-1)+NFRSTP(J-1)
      NTESTP=NADDP(NTALI)+NFRSTP(NTALI)
      NTESTP=NTESTP*NRAD
      IF (NTESTP.NE.NPLPRM) THEN
        WRITE (iunout,*) 'PARAMETER ERROR DETECTED IN SETPRM: NPLPRM'
        WRITE (iunout,*) 'NTESTP, NPLPRM ',NTESTP,NPLPRM
        CALL EXIT_OWN(1)
      ENDIF


      IF (TRCTAL) THEN
        CALL LEER(2)
        WRITE(IUNOUT,*) 'VOLUME AVERAGED TALLIES CALCULATED IN THIS RUN'
        CALL LEER(1)
        WRITE(IUNOUT,'(A6,1X,A)') 'NO.','DESCRIPTION'
        DO ITAL=1,NTALV
          IF (LIVTALV(ITAL))
     .      WRITE (IUNOUT,'(I6,1X,A72)') ITAL,TXTTAL(1,ITAL)
        END DO

        IF (.NOT.ALL(LIVTALV)) THEN
          CALL LEER(2)
          WRITE(IUNOUT,*) 'VOLUME AVERAGED TALLIES NOT CALCULATED ',
     .                    'IN THIS RUN'
          CALL LEER(1)
          WRITE(IUNOUT,'(A6,1X,A)') 'NO.','DESCRIPTION'
          DO ITAL=1,NTALV
            IF (.NOT.LIVTALV(ITAL))
     .        WRITE (IUNOUT,'(I6,1X,A72)') ITAL,TXTTAL(1,ITAL)
          END DO
        END IF

        IF (ANY(LMISTALV)) THEN
          CALL LEER(2)
          WRITE(IUNOUT,*) 'VOLUME AVERAGED TALLIES EXPLICITELY ',
     .                'SWITCHED OFF VIA INPUT FILE '
          CALL LEER(1)
          WRITE(IUNOUT,'(A6,1X,A)') 'NO.','DESCRIPTION'
          DO ITAL=1,NTALV
            IF (LMISTALV(ITAL))
     .        WRITE (IUNOUT,'(I6,1X,A72)') ITAL,TXTTAL(1,ITAL)
          END DO
        END IF
        
        CALL LEER(2)
        WRITE(IUNOUT,*)'SURFACE AVERAGED TALLIES CALCULATED IN THIS RUN'
        CALL LEER(1)
        WRITE(IUNOUT,'(A6,1X,A)') 'NO.','DESCRIPTION'
        DO ITAL=1,NTALS
          IF (LIVTALS(ITAL))
     .      WRITE (IUNOUT,'(I6,1X,A72)') ITAL,TXTTLW(1,ITAL)
        END DO

        IF (.NOT.ALL(LIVTALS)) THEN
          CALL LEER(2)
          WRITE(IUNOUT,*) 'SURFACE AVERAGED TALLIES NOT CALCULATED ',
     .                    'IN THIS RUN'
          CALL LEER(1)
          WRITE(IUNOUT,'(A6,1X,A)') 'NO.','DESCRIPTION'
          DO ITAL=1,NTALS
            IF (.NOT.LIVTALS(ITAL))
     .        WRITE (IUNOUT,'(I6,1X,A72)') ITAL,TXTTLW(1,ITAL)
          END DO
        END IF

        IF (ANY(LMISTALS)) THEN
          CALL LEER(2)
          WRITE(IUNOUT,*) 'SURFACE AVERAGED TALLIES EXPLICITELY ',
     .                'SWITCHED OFF VIA INPUT FILE '
          CALL LEER(1)
          WRITE(IUNOUT,'(A6,1X,A)') 'NO.','DESCRIPTION'
          DO ITAL=1,NTALS
            IF (LMISTALS(ITAL))
     .        WRITE (IUNOUT,'(I6,1X,A72)') ITAL,TXTTLW(1,ITAL)
          END DO
        END IF

        CALL LEER(1)
        
      END IF
C
      RETURN
      END

C ===== SOURCE: settxt.f
c  bug fix:  text(71-13) --> text(71)
c  17.3.06: txttal and txttlw added for additional tallies
C
      SUBROUTINE SETTXT

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CTEXT
      USE COUTAU

      IMPLICIT NONE

      INTEGER :: IATM, IION, IPLS, IMOL, ISPZ, IPHOT, I, J, N1,
     .           N2, N3, N4, N5, N6, N7, N8, N9, N10, N11
      CHARACTER(24) :: TEXT24
      CHARACTER(72) :: TEXT72

      TXTTAL(1,1)='PARTICLE DENSITY (ATOMS)                         '
      TXTTAL(1,2)='PARTICLE DENSITY (MOLECULES)                     '
      TXTTAL(1,3)='PARTICLE DENSITY (TEST IONS)                     '
      TXTTAL(1,4)='PARTICLE DENSITY (PHOTONS)                       '
      TXTTAL(1,5)='ENERGY DENSITY (ATOMS)                           '
      TXTTAL(1,6)='ENERGY DENSITY (MOLECULES)                       '
      TXTTAL(1,7)='ENERGY DENSITY (TEST IONS)                       '
      TXTTAL(1,8)='ENERGY DENSITY (PHOTONS)                         '
      TXTTAL(1,9)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,10)=
     . 'PARTICLE SOURCE (ATOMS) FROM ATOM-PLASMA INTERACTION        '
      TXTTAL(1,11)=
     . 'PARTICLE SOURCE (MOLECULES) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,12)=
     . 'PARTICLE SOURCE (TEST IONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,13)=
     . 'PARTICLE SOURCE (PHOTONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,14)=
     . 'PARTICLE SOURCE (BULK IONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,15)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,16)=
     . 'PARTICLE SOURCE (ATOMS) FROM MOLECULE-PLASMA INTERACTION    '
      TXTTAL(1,17)=
     . 'PARTICLE SOURCE (MOLECULES) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,18)=
     . 'PARTICLE SOURCE (TEST IONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,19)=
     . 'PARTICLE SOURCE (PHOTONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,20)=
     . 'PARTICLE SOURCE (BULK IONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,21)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,22)=
     . 'PARTICLE SOURCE (ATOMS) FROM TEST ION-PLASMA INTERACTION    '
      TXTTAL(1,23)=
     . 'PARTICLE SOURCE (MOLECULES) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,24)=
     . 'PARTICLE SOURCE (TEST IONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,25)=
     . 'PARTICLE SOURCE (PHOTONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,26)=
     . 'PARTICLE SOURCE (BULK IONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,27)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,28)=
     . 'PARTICLE SOURCE (ATOMS) FROM PHOTON-PLASMA INTERACTION      '
      TXTTAL(1,29)=
     . 'PARTICLE SOURCE (MOLECULES) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,30)=
     . 'PARTICLE SOURCE (TEST IONS) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,31)=
     . 'PARTICLE SOURCE (PHOTONS) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,32)=
     . 'PARTICLE SOURCE (BULK IONS) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,33)=
     . 'ENERGY SOURCE (ELECTRONS) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,34)=
     . 'ENERGY SOURCE (ATOMS) FROM ATOM-PLASMA INTERACTION          '
      TXTTAL(1,35)=
     . 'ENERGY SOURCE (MOLECULES) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,36)=
     . 'ENERGY SOURCE (TEST IONS) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,37)=
     . 'ENERGY SOURCE (PHOTONS) FROM ATOM-PLASMA INTERACTION        '
      TXTTAL(1,38)=
     . 'ENERGY SOURCE (BULK IONS) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,39)=
     . 'ENERGY SOURCE (ELECTRONS) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,40)=
     . 'ENERGY SOURCE (ATOMS) FROM MOLECULE-PLASMA INTERACTION      '
      TXTTAL(1,41)=
     . 'ENERGY SOURCE (MOLECULES) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,42)=
     . 'ENERGY SOURCE (TEST IONS) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,43)=
     . 'ENERGY SOURCE (PHOTONS) FROM MOLECULE-PLASMA INTERACTION    '
      TXTTAL(1,44)=
     . 'ENERGY SOURCE (BULK IONS) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,45)=
     . 'ENERGY SOURCE (ELECTRONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,46)=
     . 'ENERGY SOURCE (ATOMS) FROM TEST ION-PLASMA INTERACTION      '
      TXTTAL(1,47)=
     . 'ENERGY SOURCE (MOLECULES) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,48)=
     . 'ENERGY SOURCE (TEST IONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,49)=
     . 'ENERGY SOURCE (PHOTONS) FROM TEST ION-PLASMA INTERACTION    '
      TXTTAL(1,50)=
     . 'ENERGY SOURCE (BULK IONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,51)=
     . 'ENERGY SOURCE (ELECTRONS) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,52)=
     . 'ENERGY SOURCE (ATOMS) FROM PHOTON-PLASMA INTERACTION        '
      TXTTAL(1,53)=
     . 'ENERGY SOURCE (MOLECULES) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,54)=
     . 'ENERGY SOURCE (TEST IONS) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,55)=
     . 'ENERGY SOURCE (PHOTONS) FROM PHOTON-PLASMA INTERACTION      '
      TXTTAL(1,56)=
     . 'ENERGY SOURCE (BULK IONS) FROM PHOTON-PLASMA INTERACTION    '
C  TALLY NTALA=57 (SEE PARMMOD.F)
C        ADDIT. TRACKLENGTH ESTIMATED TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 10A
      TXTTAL(1,NTALA)=
     . 'ADDITIONAL TALLIES, TRACKLENGTH ESTIMATOR, SUBR. UPTUSR.F   '
C  TALLY NTALC=58 (SEE PARMMOD.F)
C        ADDIT. COLLISION ESTIMATED TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 10B
      TXTTAL(1,NTALC)=
     . 'ADDITIONAL TALLIES, COLLISION ESTIMATOR, SUBR. UPCUSR.F     '
C  TALLY NTALT=59 (SEE PARMMOD.F)
C        ADDIT. SNAPSHOT ESTIMATED TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 13B
      TXTTAL(1,NTALT)=
     . 'ADDITIONAL TALLIES, SNAPSHOT ESTIMATOR, SUBR. UPNUSR.F      '
C  TALLY NTALM=60 (SEE PARMMOD.F)
C        ADDIT. TALLIES FOR INTERFACING TO OTHER CODES
C        TXTTAL MAY BE OVERWRITTEN IN SUBR. INFCOP
      TXTTAL(1,NTALM)=
     . 'ADDITIONAL TALLIES FOR INTERFACING, SUBR. INFCOP.F          '
C  TALLY NTALB=61 (SEE PARMMOD.F)
C        ADDIT. TALLIES FOR ITERATIVE MODE (BGK-ITERATION)
      TXTTAL(1,NTALB)=
     . 'ADDITIONAL TALLIES FOR ITERATIVE MODE, SUBR. UPTBGK.F       '
C  TALLY NTALB=62 (SEE PARMMOD.F)
C        ADDIT. TALLIES, ALGEBRAIC EXPRESSION IN EXISTING TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 10C
      TXTTAL(1,NTALR)=
     . 'ADDITIONAL TALLIES, ALGEBRAIC EXPRESSIONS, INPUT BLOCK 10C  '
      TXTTAL(1,63)=
     . 'PARTICLE SINK (ATOMS) DUE TO GENERATION LIMIT               '
      TXTTAL(1,64)=
     . 'PARTICLE SINK (MOLECULES) DUE TO GENERATION LIMIT           '
      TXTTAL(1,65)=
     . 'PARTICLE SINK (TEST IONS) DUE TO GENERATION LIMIT           '
      TXTTAL(1,66)=
     . 'PARTICLE SINK (PHOTONS) DUE TO GENERATION LIMIT             '
      TXTTAL(1,67)=
     . 'ENERGY SINK (ATOMS) DUE TO GENERATION LIMIT                 '
      TXTTAL(1,68)=
     . 'ENERGY SINK (MOLECULES) DUE TO GENERATION LIMIT             '
      TXTTAL(1,69)=
     . 'ENERGY SINK (TEST IONS) DUE TO GENERATION LIMIT             '
      TXTTAL(1,70)=
     . 'ENERGY SINK (PHOTONS) DUE TO GENERATION LIMIT               '
      TXTTAL(1,71)=
     . 'MOMENTUM SINK (ATOMS) DUE TO GENERATION LIMIT               '
      TXTTAL(1,72)=
     . 'MOMENTUM SINK (MOLECULES) DUE TO GENERATION LIMIT           '
      TXTTAL(1,73)=
     . 'MOMENTUM SINK (TEST IONS) DUE TO GENERATION LIMIT           '
      TXTTAL(1,74)=
     . 'MOMENTUM SINK (PHOTONS) DUE TO GENERATION LIMIT             '
      TXTTAL(1,75)=
     . 'PRIMARY PARTICLE SOURCE (ATOMS) FROM PLASMA INTERACTIONS    '
      TXTTAL(1,76)=
     . 'PRIMARY PARTICLE SOURCE (MOLECULES) FROM PLASMA INTERACTIONS'
      TXTTAL(1,77)=
     . 'PRIMARY PARTICLE SOURCE (TEST IONS) FROM PLASMA INTERACTIONS'
      TXTTAL(1,78)=
     . 'PRIMARY PARTICLE SOURCE (PHOTONS) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,79)=
     . 'PRIMARY PARTICLE SOURCE (BULK IONS) FROM PLASMA INTERACTIONS'
      TXTTAL(1,80)=
     . 'PRIMARY ENERGY SOURCE (ATOMS) FROM PLASMA INTERACTIONS      '
      TXTTAL(1,81)=
     . 'PRIMARY ENERGY SOURCE (MOLECULES) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,82)=
     . 'PRIMARY ENERGY SOURCE (TEST IONS) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,83)=
     . 'PRIMARY ENERGY SOURCE (PHOTONS) FROM PLASMA INTERACTIONS    '
      TXTTAL(1,84)=
     . 'PRIMARY ENERGY SOURCE (BULK IONS) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,85)=
     . 'MOMENTUM DENSITY, X-DIRECTION (ATOMS)                       '
      TXTTAL(1,86)=
     . 'MOMENTUM DENSITY, X-DIRECTION (MOLECULES)                   '
      TXTTAL(1,87)=
     . 'MOMENTUM DENSITY, X-DIRECTION (TEST IONS)                   '
      TXTTAL(1,88)=
     . 'MOMENTUM DENSITY, X-DIRECTION (PHOTONS)                     '
      TXTTAL(1,89)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (ATOMS)                       '
      TXTTAL(1,90)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (MOLECULES)                   '
      TXTTAL(1,91)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (TEST IONS)                   '
      TXTTAL(1,92)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (PHOTONS)                     '
      TXTTAL(1,93)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (ATOMS)                       '
      TXTTAL(1,94)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (MOLECULES)                   '
      TXTTAL(1,95)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (TEST IONS)                   '
      TXTTAL(1,96)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (PHOTONS)                     '
      TXTTAL(1,97)=
     . 'MOMENTUM SOURCE (BULK IONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,98)=
     . 'MOMENTUM SOURCE (BULK IONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,99)=
     . 'MOMENTUM SOURCE (BULK IONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,100)=
     . 'MOMENTUM SOURCE (BULK IONS) FROM PHOTON-PLASMA INTERACTION  '
C
      DO 1 J=1,NTALV
        DO 1 I=2,N1MX
          TEXT72=TXTTAL(1,J)
          TXTTAL(I,J)=TEXT72
1     CONTINUE
C
      TXTUNT(1,1)='CM**-3                  '
      TXTUNT(1,2)='CM**-3                  '
      TXTUNT(1,3)='CM**-3                  '
      TXTUNT(1,4)='CM**-3                  '
      TXTUNT(1,5)='EV*CM**-3               '
      TXTUNT(1,6)='EV*CM**-3               '
      TXTUNT(1,7)='EV*CM**-3               '
      TXTUNT(1,8)='EV*CM**-3               '
      TXTUNT(1,9)='AMP*CM**-3              '
      TXTUNT(1,10)='AMP*CM**-3              '
      TXTUNT(1,11)='AMP*CM**-3              '
      TXTUNT(1,12)='AMP*CM**-3              '
      TXTUNT(1,13)='AMP*CM**-3              '
      TXTUNT(1,14)='AMP*CM**-3              '
      TXTUNT(1,15)='AMP*CM**-3              '
      TXTUNT(1,16)='AMP*CM**-3              '
      TXTUNT(1,17)='AMP*CM**-3              '
      TXTUNT(1,18)='AMP*CM**-3              '
      TXTUNT(1,19)='AMP*CM**-3              '
      TXTUNT(1,20)='AMP*CM**-3              '
      TXTUNT(1,21)='AMP*CM**-3              '
      TXTUNT(1,22)='AMP*CM**-3              '
      TXTUNT(1,23)='AMP*CM**-3              '
      TXTUNT(1,24)='AMP*CM**-3              '
      TXTUNT(1,25)='AMP*CM**-3              '
      TXTUNT(1,26)='AMP*CM**-3              '
      TXTUNT(1,27)='AMP*CM**-3              '
      TXTUNT(1,28)='AMP*CM**-3              '
      TXTUNT(1,29)='AMP*CM**-3              '
      TXTUNT(1,30)='AMP*CM**-3              '
      TXTUNT(1,31)='AMP*CM**-3              '
      TXTUNT(1,32)='AMP*CM**-3              '
      TXTUNT(1,33)='WATT*CM**-3             '
      TXTUNT(1,34)='WATT*CM**-3             '
      TXTUNT(1,35)='WATT*CM**-3             '
      TXTUNT(1,36)='WATT*CM**-3             '
      TXTUNT(1,37)='WATT*CM**-3             '
      TXTUNT(1,38)='WATT*CM**-3             '
      TXTUNT(1,39)='WATT*CM**-3             '
      TXTUNT(1,40)='WATT*CM**-3             '
      TXTUNT(1,41)='WATT*CM**-3             '
      TXTUNT(1,42)='WATT*CM**-3             '
      TXTUNT(1,43)='WATT*CM**-3             '
      TXTUNT(1,44)='WATT*CM**-3             '
      TXTUNT(1,45)='WATT*CM**-3             '
      TXTUNT(1,46)='WATT*CM**-3             '
      TXTUNT(1,47)='WATT*CM**-3             '
      TXTUNT(1,48)='WATT*CM**-3             '
      TXTUNT(1,49)='WATT*CM**-3             '
      TXTUNT(1,50)='WATT*CM**-3             '
      TXTUNT(1,51)='WATT*CM**-3             '
      TXTUNT(1,52)='WATT*CM**-3             '
      TXTUNT(1,53)='WATT*CM**-3             '
      TXTUNT(1,54)='WATT*CM**-3             '
      TXTUNT(1,55)='WATT*CM**-3             '
      TXTUNT(1,56)='WATT*CM**-3             '
      TXTUNT(1,NTALA)='TO BE READ              '
      TXTUNT(1,NTALC)='TO BE READ              '
      TXTUNT(1,NTALM)='TO BE DEFINED IN INFCOP '
      TXTUNT(1,NTALR)='TO BE READ              '
      TXTUNT(1,NTALB)='TO BE DEFINED IN BGK    '
C  GENERATION LIMIT TALLIES
      TXTUNT(1,63)='AMP*CM**-3              '
      TXTUNT(1,64)='AMP*CM**-3              '
      TXTUNT(1,65)='AMP*CM**-3              '
      TXTUNT(1,66)='AMP*CM**-3              '
      TXTUNT(1,67)='EV*CM**-3               '
      TXTUNT(1,68)='EV*CM**-3               '
      TXTUNT(1,69)='EV*CM**-3               '
      TXTUNT(1,70)='EV*CM**-3               '
      TXTUNT(1,71)='CM/S*CM**-3             '
      TXTUNT(1,72)='CM/S*CM**-3             '
      TXTUNT(1,73)='CM/S*CM**-3             '
      TXTUNT(1,74)='CM/S*CM**-3             '
      TXTUNT(1,75)='AMP*CM**-3              '
      TXTUNT(1,76)='AMP*CM**-3              '
      TXTUNT(1,77)='AMP*CM**-3              '
      TXTUNT(1,78)='AMP*CM**-3              '
      TXTUNT(1,79)='AMP*CM**-3              '
      TXTUNT(1,80)='WATT*CM**-3             '
      TXTUNT(1,81)='WATT*CM**-3             '
      TXTUNT(1,82)='WATT*CM**-3             '
      TXTUNT(1,83)='WATT*CM**-3             '
      TXTUNT(1,84)='WATT*CM**-3             '
      TXTUNT(1,85)='G*CM/SEC*CM**-3         '
      TXTUNT(1,86)='G*CM/SEC*CM**-3         '
      TXTUNT(1,87)='G*CM/SEC*CM**-3         '
      TXTUNT(1,88)='G*CM/SEC*CM**-3         '
      TXTUNT(1,89)='G*CM/SEC*CM**-3         '
      TXTUNT(1,90)='G*CM/SEC*CM**-3         '
      TXTUNT(1,91)='G*CM/SEC*CM**-3         '
      TXTUNT(1,92)='G*CM/SEC*CM**-3         '
      TXTUNT(1,93)='G*CM/SEC*CM**-3         '
      TXTUNT(1,94)='G*CM/SEC*CM**-3         '
      TXTUNT(1,95)='G*CM/SEC*CM**-3         '
      TXTUNT(1,96)='G*CM/SEC*CM**-3         '
      TXTUNT(1,97)='AMP*CM**-3              '
      TXTUNT(1,98)='AMP*CM**-3              '
      TXTUNT(1,99)='AMP*CM**-3              '
      TXTUNT(1,100)='AMP*CM**-3             '
      DO 2 J=1,NTALV
        DO 2 I=2,N1MX
          TEXT24=TXTUNT(1,J)
          TXTUNT(I,J)=TEXT24
2     CONTINUE
C
      TXTTLW(1,1)='PARTICLE FLUX, INCIDENT, ATOMS                   '
      TXTTLW(1,2)='PARTICLE FLUX, EMITTED, ATS. => ATOMS            '
      TXTTLW(1,3)='PARTICLE FLUX, EMITTED, MLS. => ATOMS            '
      TXTTLW(1,4)='PARTICLE FLUX, EMITTED, T.I. => ATOMS            '
      TXTTLW(1,5)='PARTICLE FLUX, EMITTED, PHS. => ATOMS            '
      TXTTLW(1,6)='PARTICLE FLUX, EMITTED, B.I. => ATOMS            '
      TXTTLW(1,7)='PARTICLE FLUX, INCIDENT, MOLECULES               '
      TXTTLW(1,8)='PARTICLE FLUX, EMITTED, ATS. => MOLECULES        '
      TXTTLW(1,9)='PARTICLE FLUX, EMITTED, MLS. => MOLECULES        '
      TXTTLW(1,10)='PARTICLE FLUX, EMITTED, T.I. => MOLECULES        '
      TXTTLW(1,11)='PARTICLE FLUX, EMITTED, PHS. => MOLECULES        '
      TXTTLW(1,12)='PARTICLE FLUX, EMITTED, B.I. => MOLECULES        '
      TXTTLW(1,13)='PARTICLE FLUX, INCIDENT, TEST IONS               '
      TXTTLW(1,14)='PARTICLE FLUX, EMITTED, ATS. => TEST IONS        '
      TXTTLW(1,15)='PARTICLE FLUX, EMITTED, MLS. => TEST IONS        '
      TXTTLW(1,16)='PARTICLE FLUX, EMITTED, T.I. => TEST IONS        '
      TXTTLW(1,17)='PARTICLE FLUX, EMITTED, PHS. => TEST IONS        '
      TXTTLW(1,18)='PARTICLE FLUX, EMITTED, B.I. => TEST IONS        '
      TXTTLW(1,19)='PARTICLE FLUX, INCIDENT, PHOTONS                 '
      TXTTLW(1,20)='PARTICLE FLUX, EMITTED, ATS. => PHOTONS          '
      TXTTLW(1,21)='PARTICLE FLUX, EMITTED, MLS. => PHOTONS          '
      TXTTLW(1,22)='PARTICLE FLUX, EMITTED, T.I. => PHOTONS          '
      TXTTLW(1,23)='PARTICLE FLUX, EMITTED, PHS. => PHOTONS          '
      TXTTLW(1,24)='PARTICLE FLUX, EMITTED, B.I. => PHOTONS          '
      TXTTLW(1,25)='PARTICLE FLUX, INCIDENT, BULK IONS               '

      TXTTLW(1,26)='ENERGY FLUX, INCIDENT, ATOMS                     '
      TXTTLW(1,27)='ENERGY FLUX, EMITTED, ATS. => ATOMS              '
      TXTTLW(1,28)='ENERGY FLUX, EMITTED, MLS. => ATOMS              '
      TXTTLW(1,29)='ENERGY FLUX, EMITTED, T.I. => ATOMS              '
      TXTTLW(1,30)='ENERGY FLUX, EMITTED, PHS. => ATOMS              '
      TXTTLW(1,31)='ENERGY FLUX, EMITTED, B.I. => ATOMS              '
      TXTTLW(1,32)='ENERGY FLUX, INCIDENT, MOLECULES                 '
      TXTTLW(1,33)='ENERGY FLUX, EMITTED, ATS. => MOLECULES          '
      TXTTLW(1,34)='ENERGY FLUX, EMITTED, MLS. => MOLECULES          '
      TXTTLW(1,35)='ENERGY FLUX, EMITTED, T.I. => MOLECULES          '
      TXTTLW(1,36)='ENERGY FLUX, EMITTED, PHS. => MOLECULES          '
      TXTTLW(1,37)='ENERGY FLUX, EMITTED, B.I. => MOLECULES          '
      TXTTLW(1,38)='ENERGY FLUX, INCIDENT, TEST IONS                 '
      TXTTLW(1,39)='ENERGY FLUX, EMITTED, ATS. => TEST IONS          '
      TXTTLW(1,40)='ENERGY FLUX, EMITTED, MLS. => TEST IONS          '
      TXTTLW(1,41)='ENERGY FLUX, EMITTED, T.I. => TEST IONS          '
      TXTTLW(1,42)='ENERGY FLUX, EMITTED, PHS. => TEST IONS          '
      TXTTLW(1,43)='ENERGY FLUX, EMITTED, B.I. => TEST IONS          '
      TXTTLW(1,44)='ENERGY FLUX, INCIDENT, PHOTONS                   '
      TXTTLW(1,45)='ENERGY FLUX, EMITTED, ATS. => PHOTONS            '
      TXTTLW(1,46)='ENERGY FLUX, EMITTED, MLS. => PHOTONS            '
      TXTTLW(1,47)='ENERGY FLUX, EMITTED, T.I. => PHOTONS            '
      TXTTLW(1,48)='ENERGY FLUX, EMITTED, PHS. => PHOTONS            '
      TXTTLW(1,49)='ENERGY FLUX, EMITTED, B.I. => PHOTONS            '
      TXTTLW(1,50)='ENERGY FLUX, INCIDENT, BULK IONS                 '

      TXTTLW(1,51)='SPUTTERED FLUX BY INCIDENT ATOMS                 '
      TXTTLW(1,52)='SPUTTERED FLUX BY INCIDENT MOLECULES             '
      TXTTLW(1,53)='SPUTTERED FLUX BY INCIDENT TEST IONS             '
      TXTTLW(1,54)='SPUTTERED FLUX BY INCIDENT PHOTONS               '
      TXTTLW(1,55)='SPUTTERED FLUX BY INCIDENT BULK IONS             '
      TXTTLW(1,56)='SPUTTERED FLUX, TOTAL                            '
C  TALLY NTLSA=57 (SEE PARMMOD.F)
C   ADDIT. TALLIES, SUBR. UPSUSR.F
C   TXTTLW IS OVERWRITTEN BY INPUT BLOCK 10D
      TXTTLW(1,57)='ADDITIONAL SURFACE TALLY, SUBR. UPSUSR.F         '
C  TALLY NTLSR=58 (SEE PARMMOD.F)
C   ADDIT. TALLIES, ALGEBRAIC EXPRESSION IN EXISTING TALLIES
C   TXTTLW IS OVERWRITTEN BY INPUT BLOCK 10E
      TXTTLW(1,58)='ALGEBRAIC EXPRESSION IN SURFACE AVERAGED TALLIES '
      TXTTLW(1,59)='PUMPED FLUX BY SPECIES                           '
C
      DO J=1,NTALS
        TXTTLW(2:N2MX,J)=TXTTLW(1,J)
      END DO
C
      TXTUNW(1,1)='AMP                     '
      TXTUNW(1,2)='AMP                     '
      TXTUNW(1,3)='AMP                     '
      TXTUNW(1,4)='AMP                     '
      TXTUNW(1,5)='AMP                     '
      TXTUNW(1,6)='AMP                     '
      TXTUNW(1,7)='AMP                     '
      TXTUNW(1,8)='AMP                     '
      TXTUNW(1,9)='AMP                     '
      TXTUNW(1,10)='AMP                     '
      TXTUNW(1,11)='AMP                     '
      TXTUNW(1,12)='AMP                     '
      TXTUNW(1,13)='AMP                     '
      TXTUNW(1,14)='AMP                     '
      TXTUNW(1,15)='AMP                     '
      TXTUNW(1,16)='AMP                     '
      TXTUNW(1,17)='AMP                     '
      TXTUNW(1,18)='AMP                     '
      TXTUNW(1,19)='AMP                     '
      TXTUNW(1,20)='AMP                     '
      TXTUNW(1,21)='AMP                     '
      TXTUNW(1,22)='AMP                     '
      TXTUNW(1,23)='AMP                     '
      TXTUNW(1,24)='AMP                     '
      TXTUNW(1,25)='AMP                     '

      TXTUNW(1,26)='WATT                    '
      TXTUNW(1,27)='WATT                    '
      TXTUNW(1,28)='WATT                    '
      TXTUNW(1,29)='WATT                    '
      TXTUNW(1,30)='WATT                    '
      TXTUNW(1,31)='WATT                    '
      TXTUNW(1,32)='WATT                    '
      TXTUNW(1,33)='WATT                    '
      TXTUNW(1,34)='WATT                    '
      TXTUNW(1,35)='WATT                    '
      TXTUNW(1,36)='WATT                    '
      TXTUNW(1,37)='WATT                    '
      TXTUNW(1,38)='WATT                    '
      TXTUNW(1,39)='WATT                    '
      TXTUNW(1,40)='WATT                    '
      TXTUNW(1,41)='WATT                    '
      TXTUNW(1,42)='WATT                    '
      TXTUNW(1,43)='WATT                    '
      TXTUNW(1,44)='WATT                    '
      TXTUNW(1,45)='WATT                    '
      TXTUNW(1,46)='WATT                    '
      TXTUNW(1,47)='WATT                    '
      TXTUNW(1,48)='WATT                    '
      TXTUNW(1,49)='WATT                    '
      TXTUNW(1,50)='WATT                    '

      TXTUNW(1,51)='AMP                     '
      TXTUNW(1,52)='AMP                     '
      TXTUNW(1,53)='AMP                     '
      TXTUNW(1,54)='AMP                     '
      TXTUNW(1,55)='AMP                     '
      TXTUNW(1,56)='AMP                     '
      TXTUNW(1,57)='TO BE READ              '
      TXTUNW(1,58)='TO BE READ              '
      TXTUNW(1,59)='AMP                     '
      DO J=1,NTALS
        TXTUNW(2:N2MX,J)=TXTUNW(1,J)
      END DO
C
      TXTPLS(1,1)='PLASMA TEMPERATURE                               '
      TXTPLS(1,2)='PLASMA TEMPERATURE                               '
      TXTPLS(1,3)='PLASMA DENSITY (BULK PARTICLES)                  '
      TXTPLS(1,4)='PLASMA DENSITY (BULK PARTICLES)                  '
      TXTPLS(1,5)='DRIFT VELOCITY IN X-DIRECTION (BULK IONS)        '
      TXTPLS(1,6)='DRIFT VELOCITY IN Y-DIRECTION (BULK IONS)        '
      TXTPLS(1,7)='DRIFT VELOCITY IN Z-DIRECTION (BULK IONS)        '
      TXTPLS(1,8)='MAGN. FIELD UNIT VECTOR, X DIRECTION             '
      TXTPLS(1,9)='MAGN. FIELD UNIT VECTOR, Y DIRECTION             '
      TXTPLS(1,10)='MAGN. FIELD UNIT VECTOR, Z DIRECTION             '
      TXTPLS(1,11)='MAGN. FIELD STRENGTH                             '
C     TXTPLS(1,12)='TO BE READ                                       '
      TXTPLS(1,13)='BULK ION KINETIC DRIFT ENERGY                    '
      TXTPLS(1,14)='ZONE VOLUMES                                     '
      TXTPLS(1,15)='SPACE-SPECIES WEIGHT FUNCTION                    '
      TXTPLS(1,16)='PERP. MAGN. FIELD VECTOR, X DIRECTION            '
      TXTPLS(1,17)='PERP. MAGN. FIELD VECTOR, Y DIRECTION            '
C
      DO J=1,NTALI
        IF (J.NE.12) THEN
          DO I=2,N1MX
            TEXT72=TXTPLS(1,J)
            TXTPLS(I,J)=TEXT72
          ENDDO
        ENDIF
      ENDDO
C
      TXTPUN(1,1)='EV                      '
      TXTPUN(1,2)='EV                      '
      TXTPUN(1,3)='CM**-3                  '
      TXTPUN(1,4)='CM**-3                  '
      TXTPUN(1,5)='CM/SEC                  '
      TXTPUN(1,6)='CM/SEC                  '
      TXTPUN(1,7)='CM/SEC                  '
      TXTPUN(1,8)=' ---                    '
      TXTPUN(1,9)=' ---                    '
      TXTPUN(1,10)=' ---                    '
      TXTPUN(1,11)='TESLA                   '
C     TXTPUN(1,12)='TO BE READ              '
      TXTPUN(1,13)='EV                      '
      TXTPUN(1,14)='CM**3                   '
      TXTPUN(1,15)=' ---                    '
      TXTPUN(1,16)=' ---                    '
      TXTPUN(1,17)=' ---                    '
C
      DO J=1,NTALI
        IF (J.NE.12) THEN
          DO I=2,N1MX
            TEXT24=TXTPUN(1,J)
            TXTPUN(I,J)=TEXT24
          ENDDO
        ENDIF
      ENDDO
      RETURN
C
      ENTRY STTXT1
C
C
      NFSTVI(1)=NATMI
      NFSTVI(2)=NMOLI
      NFSTVI(3)=NIONI
      NFSTVI(4)=NPHOTI
      NFSTVI(5)=NATMI
      NFSTVI(6)=NMOLI
      NFSTVI(7)=NIONI
      NFSTVI(8)=NPHOTI
      NFSTVI(9)=1
      NFSTVI(10)=NATMI
      NFSTVI(11)=NMOLI
      NFSTVI(12)=NIONI
      NFSTVI(13)=NPHOTI
      NFSTVI(14)=NPLSI
      NFSTVI(15)=1
      NFSTVI(16)=NATMI
      NFSTVI(17)=NMOLI
      NFSTVI(18)=NIONI
      NFSTVI(19)=NPHOTI
      NFSTVI(20)=NPLSI
      NFSTVI(21)=1
      NFSTVI(22)=NATMI
      NFSTVI(23)=NMOLI
      NFSTVI(24)=NIONI
      NFSTVI(25)=NPHOTI
      NFSTVI(26)=NPLSI
      NFSTVI(27)=1
      NFSTVI(28)=NATMI
      NFSTVI(29)=NMOLI
      NFSTVI(30)=NIONI
      NFSTVI(31)=NPHOTI
      NFSTVI(32)=NPLSI
      NFSTVI(33)=1
      NFSTVI(34)=1
      NFSTVI(35)=1
      NFSTVI(36)=1
      NFSTVI(37)=1
      NFSTVI(38)=1
      NFSTVI(39)=1
      NFSTVI(40)=1
      NFSTVI(41)=1
      NFSTVI(42)=1
      NFSTVI(43)=1
      NFSTVI(44)=1
      NFSTVI(45)=1
      NFSTVI(46)=1
      NFSTVI(47)=1
      NFSTVI(48)=1
      NFSTVI(49)=1
      NFSTVI(50)=1
      NFSTVI(51)=1
      NFSTVI(52)=1
      NFSTVI(53)=1
      NFSTVI(54)=1
      NFSTVI(55)=1
      NFSTVI(56)=1
      NFSTVI(NTALA)=NADVI
      NFSTVI(NTALC)=NCLVI
      NFSTVI(NTALT)=NSNVI
      NFSTVI(NTALM)=NCPVI
C     NFSTVI(NTALB) IS DEFINED IN SUBR. XSECT...
      NFSTVI(NTALB)=0
      NFSTVI(NTALR)=NALVI
      NFSTVI(75)=NATMI
      NFSTVI(76)=NMOLI
      NFSTVI(77)=NIONI
      NFSTVI(78)=NPHOTI
      NFSTVI(79)=NPLSI
      NFSTVI(80)=1
      NFSTVI(81)=1
      NFSTVI(82)=1
      NFSTVI(83)=1
      NFSTVI(84)=1
      NFSTVI(85)=NATMI
      NFSTVI(86)=NMOLI
      NFSTVI(87)=NIONI
      NFSTVI(88)=NPHOTI
      NFSTVI(89)=NATMI
      NFSTVI(90)=NMOLI
      NFSTVI(91)=NIONI
      NFSTVI(92)=NPHOTI
      NFSTVI(93)=NATMI
      NFSTVI(94)=NMOLI
      NFSTVI(95)=NIONI
      NFSTVI(96)=NPHOTI
      NFSTVI(97)=NPLSI
      NFSTVI(98)=NPLSI
      NFSTVI(99)=NPLSI
      NFSTVI(100)=NPLSI

C
C
      NFSTWI(1)=NATMI
      NFSTWI(2)=NATMI
      NFSTWI(3)=NATMI
      NFSTWI(4)=NATMI
      NFSTWI(5)=NATMI
      NFSTWI(6)=NATMI
      NFSTWI(7)=NMOLI
      NFSTWI(8)=NMOLI
      NFSTWI(9)=NMOLI
      NFSTWI(10)=NMOLI
      NFSTWI(11)=NMOLI
      NFSTWI(12)=NMOLI
      NFSTWI(13)=NIONI
      NFSTWI(14)=NIONI
      NFSTWI(15)=NIONI
      NFSTWI(16)=NIONI
      NFSTWI(17)=NIONI
      NFSTWI(18)=NIONI
      NFSTWI(19)=NPHOTI
      NFSTWI(20)=NPHOTI
      NFSTWI(21)=NPHOTI
      NFSTWI(22)=NPHOTI
      NFSTWI(23)=NPHOTI
      NFSTWI(24)=NPHOTI
      NFSTWI(25)=NPLSI

      NFSTWI(26)=NATMI
      NFSTWI(27)=NATMI
      NFSTWI(28)=NATMI
      NFSTWI(29)=NATMI
      NFSTWI(30)=NATMI
      NFSTWI(31)=NATMI
      NFSTWI(32)=NMOLI
      NFSTWI(33)=NMOLI
      NFSTWI(34)=NMOLI
      NFSTWI(35)=NMOLI
      NFSTWI(36)=NMOLI
      NFSTWI(37)=NMOLI
      NFSTWI(38)=NIONI
      NFSTWI(39)=NIONI
      NFSTWI(40)=NIONI
      NFSTWI(41)=NIONI
      NFSTWI(42)=NIONI
      NFSTWI(43)=NIONI
      NFSTWI(44)=NPHOTI
      NFSTWI(45)=NPHOTI
      NFSTWI(46)=NPHOTI
      NFSTWI(47)=NPHOTI
      NFSTWI(48)=NPHOTI
      NFSTWI(49)=NPHOTI
      NFSTWI(50)=NPLSI

      NFSTWI(51)=NATMI
      NFSTWI(52)=NMOLI
      NFSTWI(53)=NIONI
      NFSTWI(54)=NPHOTI
      NFSTWI(55)=NPLSI
      NFSTWI(56)=1
      NFSTWI(NTLSA)=NADSI
      NFSTWI(NTLSR)=NALSI
      NFSTWI(NTALS)=NSPTOT
C
C
      NFSTPI(1)=1
      NFSTPI(2)=NPLSI
      NFSTPI(3)=1
      NFSTPI(4)=NPLSI
      NFSTPI(5)=NPLSI
      NFSTPI(6)=NPLSI
      NFSTPI(7)=NPLSI
      NFSTPI(8)=1
      NFSTPI(9)=1
      NFSTPI(10)=1
      NFSTPI(11)=1
      NFSTPI(12)=NAINI
      NFSTPI(13)=NPLSI
      NFSTPI(14)=1
      NFSTPI(15)=NATMI+NMOLI+NIONI
      NFSTPI(16)=1
      NFSTPI(17)=1
C
C  INITIALISE SPECIES ARRAYS FOR VOLUME TALLIES

      N1=NPHOTI
      N2=N1+NATMI
      N3=N2+NMOLI
      N4=N3+NIONI
      N5=N4+NPLSI
      N6=N5+NADVI
      N7=N6+NALVI
      N8=N7+NCLVI
      N9=N8+NCPVI
      N10=N9+NBGVI
      N11=N10+NSNVI

      NSPAN(1)=N1+1
      NSPAN(2)=N2+1
      NSPAN(3)=N3+1
      NSPAN(4)=1
      NSPAN(5)=N1+1
      NSPAN(6)=N2+1
      NSPAN(7)=N3+1
      NSPAN(8)=1
      NSPAN(9)=0
      NSPAN(10)=N1+1
      NSPAN(11)=N2+1
      NSPAN(12)=N3+1
      NSPAN(13)=1
      NSPAN(14)=N4+1
      NSPAN(15)=0
      NSPAN(16)=N1+1
      NSPAN(17)=N2+1
      NSPAN(18)=N3+1
      NSPAN(19)=1
      NSPAN(20)=N4+1
      NSPAN(21)=0
      NSPAN(22)=N1+1
      NSPAN(23)=N2+1
      NSPAN(24)=N2+1
      NSPAN(25)=1
      NSPAN(26)=N4+1
      NSPAN(27)=0
      NSPAN(28)=N1+1
      NSPAN(29)=N2+1
      NSPAN(30)=N2+1
      NSPAN(31)=1
      NSPAN(32)=N4+1
      NSPAN(33)=0
      NSPAN(34)=0
      NSPAN(35)=0
      NSPAN(36)=0
      NSPAN(37)=0
      NSPAN(38)=0
      NSPAN(39)=0
      NSPAN(40)=0
      NSPAN(41)=0
      NSPAN(42)=0
      NSPAN(43)=0
      NSPAN(44)=0
      NSPAN(45)=0
      NSPAN(46)=0
      NSPAN(47)=0
      NSPAN(48)=0
      NSPAN(49)=0
      NSPAN(50)=0
      NSPAN(51)=0
      NSPAN(52)=0
      NSPAN(53)=0
      NSPAN(54)=0
      NSPAN(55)=0
      NSPAN(56)=0
      NSPAN(NTALA)=N5+1
      NSPAN(NTALC)=N7+1
      NSPAN(NTALT)=N10+1
      NSPAN(NTALM)=N8+1
      NSPAN(NTALB)=N9+1
      NSPAN(NTALR)=N6+1
C  GENERATION LIMIT TALLIES
      NSPAN(63)=N1+1
      NSPAN(64)=N2+1
      NSPAN(65)=N3+1
      NSPAN(66)=1
      NSPAN(67)=N1+1
      NSPAN(68)=N2+1
      NSPAN(69)=N3+1
      NSPAN(70)=1
      NSPAN(71)=N1+1
      NSPAN(72)=N2+1
      NSPAN(73)=N3+1
      NSPAN(74)=1
      NSPAN(75)=N1+1
      NSPAN(76)=N2+1
      NSPAN(77)=N3+1
      NSPAN(78)=1
      NSPAN(79)=N4+1
      NSPAN(80)=0
      NSPAN(81)=0
      NSPAN(82)=0
      NSPAN(83)=0
      NSPAN(84)=0
      NSPAN(85)=N1+1
      NSPAN(86)=N2+1
      NSPAN(87)=N3+1
      NSPAN(88)=1
      NSPAN(89)=N1+1
      NSPAN(90)=N2+1
      NSPAN(91)=N3+1
      NSPAN(92)=1
      NSPAN(93)=N1+1
      NSPAN(94)=N2+1
      NSPAN(95)=N3+1
      NSPAN(96)=1
      NSPAN(97)=N4+1
      NSPAN(98)=N4+1
      NSPAN(99)=N4+1
      NSPAN(100)=N4+1

      NSPEN(1)=N2
      NSPEN(2)=N3
      NSPEN(3)=N4
      NSPEN(4)=N1
      NSPEN(5)=N2
      NSPEN(6)=N3
      NSPEN(7)=N4
      NSPEN(8)=N1
      NSPEN(9)=0
      NSPEN(10)=N2
      NSPEN(11)=N3
      NSPEN(12)=N4
      NSPEN(13)=N1
      NSPEN(14)=N5
      NSPEN(15)=0
      NSPEN(16)=N2
      NSPEN(17)=N3
      NSPEN(18)=N4
      NSPEN(19)=N1
      NSPEN(20)=N5
      NSPEN(21)=0
      NSPEN(22)=N2
      NSPEN(23)=N3
      NSPEN(24)=N4
      NSPEN(25)=N1
      NSPEN(26)=N5
      NSPEN(27)=0
      NSPEN(28)=N2
      NSPEN(29)=N3
      NSPEN(30)=N4
      NSPEN(31)=N1
      NSPEN(32)=N5
      NSPEN(33)=0
      NSPEN(34)=0
      NSPEN(35)=0
      NSPEN(36)=0
      NSPEN(37)=0
      NSPEN(38)=0
      NSPEN(39)=0
      NSPEN(40)=0
      NSPEN(41)=0
      NSPEN(42)=0
      NSPEN(43)=0
      NSPEN(44)=0
      NSPEN(45)=0
      NSPEN(46)=0
      NSPEN(47)=0
      NSPEN(48)=0
      NSPEN(49)=0
      NSPEN(50)=0
      NSPEN(51)=0
      NSPEN(52)=0
      NSPEN(53)=0
      NSPEN(54)=0
      NSPEN(55)=0
      NSPEN(56)=0
      NSPEN(NTALA)=N6
      NSPEN(NTALC)=N8
      NSPEN(NTALT)=N11
      NSPEN(NTALM)=N9
      NSPEN(NTALB)=N10
      NSPEN(NTALR)=N7
C  GENERATION LIMIT TALLIES
      NSPEN(63)=N2
      NSPEN(64)=N3
      NSPEN(65)=N4
      NSPEN(66)=N1
      NSPEN(67)=N2
      NSPEN(68)=N3
      NSPEN(69)=N4
      NSPEN(70)=N1
      NSPEN(71)=N2
      NSPEN(72)=N3
      NSPEN(73)=N4
      NSPEN(74)=N1
      NSPEN(75)=N2
      NSPEN(76)=N3
      NSPEN(77)=N4
      NSPEN(78)=N1
      NSPEN(79)=N5
      NSPEN(80)=0
      NSPEN(81)=0
      NSPEN(82)=0
      NSPEN(83)=0
      NSPEN(84)=0
      NSPEN(85)=N2
      NSPEN(86)=N3
      NSPEN(87)=N4
      NSPEN(88)=N1
      NSPEN(89)=N2
      NSPEN(90)=N3
      NSPEN(91)=N4
      NSPEN(92)=N1
      NSPEN(93)=N2
      NSPEN(94)=N3
      NSPEN(95)=N4
      NSPEN(96)=N1
      NSPEN(97)=N5
      NSPEN(98)=N5
      NSPEN(99)=N5
      NSPEN(100)=N5

      DO IPHOT=1,NPHOTI
        ISPZ=IPHOT
        TXTSPC(IPHOT,4)=TEXTS(ISPZ)
        TXTSPC(IPHOT,8)=TEXTS(ISPZ)
        TXTSPC(IPHOT,13)=TEXTS(ISPZ)
        TXTSPC(IPHOT,19)=TEXTS(ISPZ)
        TXTSPC(IPHOT,25)=TEXTS(ISPZ)
        TXTSPC(IPHOT,31)=TEXTS(ISPZ)
        TXTSPC(IPHOT,66)=TEXTS(ISPZ)
        TXTSPC(IPHOT,70)=TEXTS(ISPZ)
        TXTSPC(IPHOT,74)=TEXTS(ISPZ)
        TXTSPC(IPHOT,78)=TEXTS(ISPZ)
        TXTSPC(IPHOT,88)=TEXTS(ISPZ)
        TXTSPC(IPHOT,92)=TEXTS(ISPZ)
        TXTSPC(IPHOT,96)=TEXTS(ISPZ)
        TXTSPC(IPHOT,100)=TEXTS(ISPZ)
      END DO

      DO 10 IATM=1,NATMI
        ISPZ=NSPH+IATM
        TXTSPC(IATM,1)=TEXTS(ISPZ)
        TXTSPC(IATM,5)=TEXTS(ISPZ)
        TXTSPC(IATM,10)=TEXTS(ISPZ)
        TXTSPC(IATM,16)=TEXTS(ISPZ)
        TXTSPC(IATM,22)=TEXTS(ISPZ)
        TXTSPC(IATM,28)=TEXTS(ISPZ)
        TXTSPC(IATM,63)=TEXTS(ISPZ)
        TXTSPC(IATM,67)=TEXTS(ISPZ)
        TXTSPC(IATM,71)=TEXTS(ISPZ)
        TXTSPC(IATM,75)=TEXTS(ISPZ)
        TXTSPC(IATM,85)=TEXTS(ISPZ)
        TXTSPC(IATM,89)=TEXTS(ISPZ)
        TXTSPC(IATM,93)=TEXTS(ISPZ)
        TXTSPC(IATM,97)=TEXTS(ISPZ)
10    CONTINUE
C
      DO 20 IMOL=1,NMOLI
        ISPZ=NSPA+IMOL
        TXTSPC(IMOL,2)=TEXTS(ISPZ)
        TXTSPC(IMOL,6)=TEXTS(ISPZ)
        TXTSPC(IMOL,11)=TEXTS(ISPZ)
        TXTSPC(IMOL,17)=TEXTS(ISPZ)
        TXTSPC(IMOL,23)=TEXTS(ISPZ)
        TXTSPC(IMOL,29)=TEXTS(ISPZ)
        TXTSPC(IMOL,64)=TEXTS(ISPZ)
        TXTSPC(IMOL,68)=TEXTS(ISPZ)
        TXTSPC(IMOL,72)=TEXTS(ISPZ)
        TXTSPC(IMOL,76)=TEXTS(ISPZ)
        TXTSPC(IMOL,86)=TEXTS(ISPZ)
        TXTSPC(IMOL,90)=TEXTS(ISPZ)
        TXTSPC(IMOL,94)=TEXTS(ISPZ)
        TXTSPC(IMOL,98)=TEXTS(ISPZ)
20    CONTINUE
C
      DO 30 IION=1,NIONI
        ISPZ=NSPAM+IION
        TXTSPC(IION,3)=TEXTS(ISPZ)
        TXTSPC(IION,7)=TEXTS(ISPZ)
        TXTSPC(IION,12)=TEXTS(ISPZ)
        TXTSPC(IION,18)=TEXTS(ISPZ)
        TXTSPC(IION,24)=TEXTS(ISPZ)
        TXTSPC(IION,30)=TEXTS(ISPZ)
        TXTSPC(IION,65)=TEXTS(ISPZ)
        TXTSPC(IION,69)=TEXTS(ISPZ)
        TXTSPC(IION,73)=TEXTS(ISPZ)
        TXTSPC(IION,77)=TEXTS(ISPZ)
        TXTSPC(IION,87)=TEXTS(ISPZ)
        TXTSPC(IION,91)=TEXTS(ISPZ)
        TXTSPC(IION,95)=TEXTS(ISPZ)
        TXTSPC(IION,99)=TEXTS(ISPZ)
30    CONTINUE
C
      DO 40 IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTSPC(IPLS,14)=TEXTS(ISPZ)
        TXTSPC(IPLS,20)=TEXTS(ISPZ)
        TXTSPC(IPLS,26)=TEXTS(ISPZ)
        TXTSPC(IPLS,32)=TEXTS(ISPZ)
        TXTSPC(IPLS,79)=TEXTS(ISPZ)
40    CONTINUE
C
      TXTSPC(1,9)='ELECTRONS               '
      TXTSPC(1,15)='ELECTRONS               '
      TXTSPC(1,21)='ELECTRONS               '
      TXTSPC(1,27)='ELECTRONS               '
      TXTSPC(1,33)='ELECTRONS               '
      TXTSPC(1,39)='ELECTRONS               '
      TXTSPC(1,45)='ELECTRONS               '
      TXTSPC(1,51)='ELECTRONS               '
C
      TXTSPC(1,34)='ATOMS                   '
      TXTSPC(1,40)='ATOMS                   '
      TXTSPC(1,46)='ATOMS                   '
      TXTSPC(1,52)='ATOMS                   '
      TXTSPC(1,80)='ATOMS                   '
C
      TXTSPC(1,35)='MOLECULES               '
      TXTSPC(1,41)='MOLECULES               '
      TXTSPC(1,47)='MOLECULES               '
      TXTSPC(1,53)='MOLECULES               '
      TXTSPC(1,81)='MOLECULES               '
C
      TXTSPC(1,36)='TEST IONS               '
      TXTSPC(1,42)='TEST IONS               '
      TXTSPC(1,48)='TEST IONS               '
      TXTSPC(1,54)='TEST IONS               '
      TXTSPC(1,82)='TEST IONS               '
C
      TXTSPC(1,37)='PHOTONS                 '
      TXTSPC(1,43)='PHOTONS                 '
      TXTSPC(1,49)='PHOTONS                 '
      TXTSPC(1,55)='PHOTONS                 '
      TXTSPC(1,83)='PHOTONS                 '
C
      TXTSPC(1,38)='BULK IONS               '
      TXTSPC(1,44)='BULK IONS               '
      TXTSPC(1,50)='BULK IONS               '
      TXTSPC(1,56)='BULK IONS               '
      TXTSPC(1,84)='BULK IONS               '
C

      DO IPHOT=1,NPHOTI
        ISPZ=IPHOT
        TXTSPW(IPHOT,19)=TEXTS(ISPZ)
        TXTSPW(IPHOT,20)=TEXTS(ISPZ)
        TXTSPW(IPHOT,21)=TEXTS(ISPZ)
        TXTSPW(IPHOT,22)=TEXTS(ISPZ)
        TXTSPW(IPHOT,23)=TEXTS(ISPZ)
        TXTSPW(IPHOT,24)=TEXTS(ISPZ)
        TXTSPW(IPHOT,44)=TEXTS(ISPZ)
        TXTSPW(IPHOT,45)=TEXTS(ISPZ)
        TXTSPW(IPHOT,46)=TEXTS(ISPZ)
        TXTSPW(IPHOT,47)=TEXTS(ISPZ)
        TXTSPW(IPHOT,48)=TEXTS(ISPZ)
        TXTSPW(IPHOT,49)=TEXTS(ISPZ)
        TXTSPW(IPHOT,54)=TEXTS(ISPZ)
      END DO

      DO IATM=1,NATMI
        ISPZ=NSPH+IATM
        TXTSPW(IATM,1)=TEXTS(ISPZ)
        TXTSPW(IATM,2)=TEXTS(ISPZ)
        TXTSPW(IATM,3)=TEXTS(ISPZ)
        TXTSPW(IATM,4)=TEXTS(ISPZ)
        TXTSPW(IATM,5)=TEXTS(ISPZ)
        TXTSPW(IATM,6)=TEXTS(ISPZ)
        TXTSPW(IATM,26)=TEXTS(ISPZ)
        TXTSPW(IATM,27)=TEXTS(ISPZ)
        TXTSPW(IATM,28)=TEXTS(ISPZ)
        TXTSPW(IATM,29)=TEXTS(ISPZ)
        TXTSPW(IATM,30)=TEXTS(ISPZ)
        TXTSPW(IATM,31)=TEXTS(ISPZ)
        TXTSPW(IATM,51)=TEXTS(ISPZ)
      END DO
C
      DO IMOL=1,NMOLI
        ISPZ=NSPA+IMOL
        TXTSPW(IMOL,7)=TEXTS(ISPZ)
        TXTSPW(IMOL,8)=TEXTS(ISPZ)
        TXTSPW(IMOL,9)=TEXTS(ISPZ)
        TXTSPW(IMOL,10)=TEXTS(ISPZ)
        TXTSPW(IMOL,11)=TEXTS(ISPZ)
        TXTSPW(IMOL,12)=TEXTS(ISPZ)
        TXTSPW(IMOL,32)=TEXTS(ISPZ)
        TXTSPW(IMOL,33)=TEXTS(ISPZ)
        TXTSPW(IMOL,34)=TEXTS(ISPZ)
        TXTSPW(IMOL,35)=TEXTS(ISPZ)
        TXTSPW(IMOL,36)=TEXTS(ISPZ)
        TXTSPW(IMOL,37)=TEXTS(ISPZ)
        TXTSPW(IMOL,52)=TEXTS(ISPZ)
      END DO
C
      DO IION=1,NIONI
        ISPZ=NSPAM+IION
        TXTSPW(IION,13)=TEXTS(ISPZ)
        TXTSPW(IION,14)=TEXTS(ISPZ)
        TXTSPW(IION,15)=TEXTS(ISPZ)
        TXTSPW(IION,16)=TEXTS(ISPZ)
        TXTSPW(IION,17)=TEXTS(ISPZ)
        TXTSPW(IION,18)=TEXTS(ISPZ)
        TXTSPW(IION,38)=TEXTS(ISPZ)
        TXTSPW(IION,39)=TEXTS(ISPZ)
        TXTSPW(IION,40)=TEXTS(ISPZ)
        TXTSPW(IION,41)=TEXTS(ISPZ)
        TXTSPW(IION,42)=TEXTS(ISPZ)
        TXTSPW(IION,43)=TEXTS(ISPZ)
        TXTSPW(IION,53)=TEXTS(ISPZ)
      END DO
C
      DO IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTSPW(IPLS,25)=TEXTS(ISPZ)
        TXTSPW(IPLS,50)=TEXTS(ISPZ)
        TXTSPW(IPLS,55)=TEXTS(ISPZ)
      END DO
C
      TXTSPW(1,42)='                        '
      TXTSPW(1,43)='                        '
      TXTSPW(1,44)='                        '
      TXTSPW(1,45)='                        '
C
      TXTPSP(1,1)='ELECTRONS               '
      TXTPSP(1,3)='ELECTRONS               '
      TXTPSP(1,8)=' ---                    '
      TXTPSP(1,9)=' ---                    '
      TXTPSP(1,10)=' ---                    '
      TXTPSP(1,11)=' ---                    '
      TXTPSP(1,14)=' ---                    '
      TXTPSP(1,16)=' ---                    '
      TXTPSP(1,17)=' ---                    '
C
C     TXTPSP(IAIN,12)='TO BE READ            '
C
      DO 50 ISPZ=1,NSPAMI
50      TXTPSP(ISPZ,15)=TEXTS(ISPZ)
C
      DO 80 IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTPSP(IPLS,2)=TEXTS(ISPZ)
        TXTPSP(IPLS,4)=TEXTS(ISPZ)
        TXTPSP(IPLS,5)=TEXTS(ISPZ)
        TXTPSP(IPLS,6)=TEXTS(ISPZ)
        TXTPSP(IPLS,7)=TEXTS(ISPZ)
80      TXTPSP(IPLS,13)=TEXTS(ISPZ)
C
      RETURN
      END







C ===== SOURCE: setup_chord_spectra.f
c  14.5.06:  bug fix: 1 line added: if nchtal.ne.1 and. nchtal.ne.3:  cycle

      subroutine setup_chord_spectra

      use precision
      use parmmod
      use cestim
      use comsig
      use comprt
      use cupd
      use ccona

      implicit none
      
      real(dp) :: c1(3), c2(3), PSIG(0:NSPZ+10)
      real(dp) :: ze, timax
      integer :: ichori, ifirst, ichrd, ipvot, nbc2, nac2, iplots, 
     .           ispc, ntot_cell, ntotsp
      
      type(spect_array), allocatable :: svestiml(:), svsmestl(:)
      TYPE(EIRENE_SPECTRUM), POINTER :: ESPEC, SSPEC
      TYPE(CELL_INFO), POINTER :: FIRST, CUR

!  FIND CELLS INTERSECTED BY CHORDS

      ze = 1._dp

      ifirst = -1
      ntot_cell = 0

      do ichori = 1, nchori

        IF (.NOT.NLSTCHR(ICHORI)) CYCLE
        IF ((NCHTAL(ICHORI) /= 1) .AND. (NCHTAL(ICHORI) /= 3)) CYCLE

        IPVOT=IPIVOT(ICHORI)
        C1(1)=XPIVOT(ICHORI)
        C1(2)=YPIVOT(ICHORI)
        C1(3)=ZPIVOT(ICHORI)
C     
        ICHRD=ICHORD(ICHORI)
        C2(1)=XCHORD(ICHORI)
        C2(2)=YCHORD(ICHORI)
        C2(3)=ZCHORD(ICHORI)
C     
        NBC2=NSPBLC(ICHORI)
        NAC2=NSPADD(ICHORI)

        ALLOCATE(TRAJ(ICHORI)%TRJ)
        TRAJ(ICHORI)%TRJ%P1 = C1
        TRAJ(ICHORI)%TRJ%P2 = C2
        TRAJ(ICHORI)%TRJ%NCOU_CELL = 0
        NULLIFY(TRAJ(ICHORI)%TRJ%CELLS)
        
        CALL LININT (IFIRST,ICHORI,C1,C2,ICHRD,IPVOT,NBC2,NAC2,ZE,
     .               PSIG,TIMAX,1,1,1,IABS(NCHENI))

        ntot_cell = ntot_cell + traj(ichori)%trj%ncou_cell
         
      end do

      IF (NTOT_CELL == 0) RETURN

!  SAVE SPECTRA SPECIFIED VIA INPUT

      IF (NADSPC > 0) THEN

        ALLOCATE(SVESTIML(NADSPC))

        DO ISPC = 1, NADSPC
          SVESTIML(ISPC)%PSPC => ESTIML(ISPC)%PSPC
        END DO
        
        DEALLOCATE(ESTIML)

        IF (ALLOCATED(SMESTL)) THEN
          ALLOCATE(SVSMESTL(NADSPC))
          DO ISPC = 1, NADSPC
            SVSMESTL(ISPC)%PSPC => SMESTL(ISPC)%PSPC
          END DO
          DEALLOCATE(SMESTL)
        END IF        

      END IF

!  set up new arrays for spectra

      NTOTSP = NADSPC + NTOT_CELL

      ALLOCATE(ESTIML(NTOTSP))
      DO ISPC = 1, NADSPC
        ESTIML(ISPC)%PSPC => SVESTIML(ISPC)%PSPC
      END DO

      IF (ALLOCATED(SVSMESTL).or.NSMSTRA.GT.0) THEN
        ALLOCATE(SMESTL(NTOTSP))
        DO ISPC = 1, NADSPC
          SMESTL(ISPC)%PSPC => SVSMESTL(ISPC)%PSPC
        END DO
      END IF
      
!  add spectra for cells along chords

      ispc = nadspc
      do ichori = 1,nchori

        IF ((NCHTAL(ICHORI) /= 1) .AND. (NCHTAL(ICHORI) /= 3)) CYCLE
         if (.not.associated(traj(ichori)%trj%cells)) cycle
         first => traj(ichori)%trj%cells
         cur => first
         do 
           allocate(espec)
           espec%isrfcll = 2
           espec%ispcsrf = cur%no_cell
           if (nchtal(ichori) == 1) then
             espec%iprtyp = 1
           else if (nchtal(ichori) == 3) then
             espec%iprtyp = 0
           else
             write (iunout,*) ' wrong chord type for spectrum '
             write (iunout,*) ' no spectrum set up for cell',cur%no_cell
             cur => cur%nextc
!pb associated with two arguments tests if both arguments point to the same target
             if (associated(cur,first)) exit
             cycle
           end if
           espec%iprsp = nspspz(ichori)
           espec%ispctyp = 1
           espec%nspc = abs(ncheni)
           espec%imetsp = 0
           espec%idirec = 1
           if (ncheni > 0) then
             espec%spcmin = emin1(ichori)
             espec%spcmax = emax1(ichori)
           else
             espec%spcmin = log10(emin1(ichori))
             espec%spcmax = log10(emax1(ichori))
           end if
           espec%esp_00 = 0._dp
           espec%spc_xplt = 0._dp
           espec%spc_yplt = 0._dp
           espec%spc_same = 0._dp
           espec%spcvx = traj(ichori)%trj%vx
           espec%spcvy = traj(ichori)%trj%vy
           espec%spcvz = traj(ichori)%trj%vz
           espec%esp_min =1.e30_dp
           espec%esp_max = -1.e30_dp
           espec%spcdel=(espec%spcmax-espec%spcmin)/real(espec%nspc,dp)
           espec%spcdeli = 1._dp / (espec%spcdel+eps60)
           allocate(espec%spc(0:espec%nspc+1))
           allocate(espec%sdv(0:espec%nspc+1))
           allocate(espec%sgm(0:espec%nspc+1))
           espec%spc(0:espec%nspc+1) = 0
           
           ispc = ispc + 1
           estiml(ispc)%pspc => espec

           if (allocated(smestl)) then
             allocate(sspec)
             allocate(sspec%spc(0:espec%nspc+1))
             allocate(sspec%sdv(0:espec%nspc+1))
             allocate(sspec%sgm(0:espec%nspc+1))
             sspec = espec
             smestl(ispc)%pspc => sspec
           end if

           cur => cur%nextc
           if (associated(cur,first)) exit
         end do
      end do

      NADSPC = ISPC

      end subroutine setup_chord_spectra
C ===== SOURCE: spltrr.f
C  AUG. 05: NCELL UPDATED FOR SPLITTING AND THEN RESET
C
      SUBROUTINE SPLTRR(IDIM,MS,NINC,*,*)
C
C  SPLITTING AND RUSSIAN ROULETTE SURFACE
C
C  IDIM  =1: RADIAL SURFACE
C        =2: POLOIDAL SURFACE
C        =3: TOROIDAL SURFACE
C        =4: ADDITIONAL SURFACE
C        =0: NOT ON ANY SURFACE
C  NINC  >0: RUSSIAN ROULETTE
C  NINC  <0: SPLITTING
C  RETURN 1: SPLIT AND CONTINUE FLIGHT
C  RETURN 2: STOP FLIGHT, BECAUSE OF RUSSIAN ROULETTE
C
C  SPLITTING AND RUSSIAN ROULETTE SURFACE MS
C
      USE PRECISION
      USE PARMMOD
      USE COMPRT
      USE COMSPL
      USE CGRID

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IDIM, MS, NINC
      REAL(DP) :: XNU, XM, YM, ZM, ZNU1, ZEP1, ZNU, X0S, Y0S, Z0S,
     .            TIMES
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: IPOLGS, MASRFS, MSURFS, J, IADD, NU, IG, NCELLS

      IF (IDIM.EQ.1) THEN
        IADD=0
      ELSEIF (IDIM.EQ.2) THEN
        IADD=N1ST
      ELSEIF (IDIM.EQ.3) THEN
        IADD=N1ST+N2ND
      ELSEIF (IDIM.EQ.4) THEN
        IADD=N1ST+N2ND+N3RD
      ELSE
      ENDIF
C
      ZNU=ABS(RNUMB(IADD+MS))
      NU=ZNU
      IG=SIGN(1._DP,RNUMB(IADD+MS))
C
C  RUSSIAN ROULETTE?
      IF(NINC*IG.GT.0) GO TO 340
C
C  NO - SPLIT PARTICLE
      X0S=X0
      Y0S=Y0
      Z0S=Z0
      TIMES=TIME
      NCELLS=NCELL
C
      X0=X0+VELX*ZT
      Y0=Y0+VELY*ZT
      Z0=Z0+VELZ*ZT
      TIME=TIME+ZT/VEL
      NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
C
C SPECIFIC FOR SPLITTING AT RADIAL SURFACES:
C
      MASRFS=MASURF
      MSURFS=MSURF
      IPOLGS=IPOLG
C
      IPOLG=IPOLGN
      MASURF=0
      MSURF=0
C
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,9)
      IF (NLSTOR) CALL STORE(200)
C
C
      NLEVEL=NLEVEL+1
C  IS SPLITTING PARAMETER AN INTEGER?
      IF (NU.NE.ZNU) THEN
        XNU=ZNU-DBLE(NU)
        IF (RANF_EIRENE( ).LT.XNU) NU=NU+1
      ENDIF
C
C   SPLIT INTO SEVERAL PARTICLES WITH REDUCED WEIGHT
C   SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
      WEIGHT=WEIGHT/ZNU
      DO 333 J=1,NPARTC
        RSPLST(NLEVEL,J)=RPST(J)
333   CONTINUE
      DO 334 J=1,MPARTC
        ISPLST(NLEVEL,J)=IPST(J)
334   CONTINUE
C  NUMBER OF NODES AT THIS LEVEL
      NODES(NLEVEL)=NU
C  RESTORE SOME PARTICLE CO-ORDINATES
      X0=X0S
      Y0=Y0S
      Z0=Z0S
      TIME=TIMES
      NCELL=NCELLS
C
C  SPECIFIC FOR SPLITTING ON RADIAL SURFACES:
      MASURF=MASRFS
      MSURF=MSURFS
      IPOLG=IPOLGS
C
C  CONTINUE THE OLD TRACK WITH REDUCED WEIGHT
      RETURN 1
C
C   RUSSIAN ROULETTE
C
340   CONTINUE
C  PROBABLITY OF DEATH
      ZNU1=1.0/ZNU
      ZEP1=RANF_EIRENE( )
C  IS THIS PARTICLE TO BE KILLED
      IF (ZEP1.LT.ZNU1) THEN
C  NO   INCREASE WEIGHT AND CONTINUE TRACKING
        WEIGHT=WEIGHT*ZNU
        RETURN 1
      ENDIF
C
C  YES   KILL PARTICLE AND STOP FLIGHT
C
      IF (NLTRC) THEN
        XM=X0+VELX*ZT
        YM=Y0+VELY*ZT
        ZM=Z0+VELZ*ZT
        CALL CHCTRC(XM,YM,ZM,16,10)
      ENDIF
      IF (NLSTOR) CALL STORE(0)
      LGPART=.FALSE.
      RETURN 2
      END
C ===== SOURCE: volume.f
C
      SUBROUTINE VOLUME (IND)

C  CALCULATE VOLUME-ELEMENTS FOR VOLUME AVERAGED TALLIES
C  THE CELL VOLUMES VOL MUST BE THOSE SEEN BY THE TESTPARTICLES
C  I.E. NOT NECESSARLY THE TRUE ONES.
C  ONE COMMON FACTOR (LENGTH OF THE CELL IN IGNORABLE
C  DIMENSION) ACTS LIKE A SCALING FACTOR FOR THESE TALLIES.
C
C  IN CASE OF THE NLTRA OPTION: IF (NLTOR):
C                               VOL = TAN(ALPHA)*XCOM*AREA*2.
C                               VOL = VOLUME OF ONE OF THE NTTRAM
C                                     CYLINDRICAL SEGMENTS
C                               AREA= AREA OF THE CELL
C                               ALPHA = 0.5*(2*PI / NTTRAM)
C                               XCOM= X - CENTER OF MASS OF AREA
C                               (I.E., NTTRAM=PI/ALPHA)
C                               IF (.NOT.NLTOR):
C                               VOL= NTTRAM TIMES THE VOLUME GIVEN
C                                    GIVEN ABOVE, IE:
C                               VOL= TAN(ALPHA)/ALPHA*2*PI*XCOM*AREA
C  IN CASE OF THE NLTRT OPTION: VOL = 2*PI*XCOM*AREA
C  (REGARDLESS OF NLTRZ,...)    VOL = VOLUME OF THE CELL, NO TOROIDAL
C                                     RESOLUTION
C                               AREA= AREA OF THE CELL
C                               XCOM= X - CENTER OF MASS OF AREA
C  NOTE: NLTRT=TRUE INTRODUCES INCONCISTENCY, SINCE PARTICLES
C        MAY SEE CYLINDER, BUT VOLUME IS COMPUTED FOR TORUS
C
C  NOTE: XCOM IS CENTER OF MASS IN TORUS SYSTEM, I.E.,
C        THE LARGE RADIUS RMTOR MUST BE ADDED TO THE
C        XCOM EVALUATED IN LOCAL SYSTEMS
C
      USE PRECISION
      USE PARMMOD
      USE COMPRT, ONLY: IUNOUT
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTETRA
      USE CTRIG

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IND
      REAL(DP), ALLOCATABLE, SAVE :: AREAP(:,:)
c slmod begin - debug
c...  /intel/Compiler/11.0/083/lib/intel64 under Linux distribution 
c     CentOS release 5.5 (Final) when running very large tetrahedron
c     grids of about 4 million objects:
c      REAL(DP), ALLOCATABLE :: AREA1(:)
c
      REAL(DP) :: AREA1(0:N1ST)
c slmod end
      REAL(DP) :: PC1(3), PC2(3), PC3(3), PC4(3)
      REAL(DP) :: AREAR, VOLSR, CAL_VOL, TWOTHIRD, VSAVE, FAC2, FAC3,
     .          PI2AT, AELL, DONE, DNULL, SY, X1, X2, Y1, XNULL, SX,
     .          ARTRIA, AR, XC, Y2, X3, Y3, X4, Y4
      INTEGER :: NCELL1, I, KP, IC4, IC, ITET, IC1, IC2, IC3, NCELLK,
     .           IT, NCELLJ, NCELL, IPP, I1ST, IRAD, IR, IFLAG, J,
     .           IP, JP, IN, K, IRP
!pb      SAVE

C     IND=1: 1-ST GRID, RAD. RESOLUTION
C     IND=2: 2-ND GRID, POL. RESOLUTION
C     IND=3: 3-RD GRID, TOR. RESOLUTION
C     IND=4: ADDITIONAL CELL REGION
C
      GOTO(100,200,300,400),IND
C
100   CONTINUE
C
      IF (.NOT.ALLOCATED(AREAP)) ALLOCATE (AREAP(N1STS,N2NDPLG))
c slmod begin - debug
c      IF (.NOT.ALLOCATED(AREA1)) ALLOCATE (AREA1(0:N1ST))
c slmod end

      DO 101 IRAD=1,NRAD
        VOL(IRAD)=0.
101   CONTINUE
      AREA1(0)=0.D0
      DO 102 I1ST=1,N1ST
        AREA1(I1ST)=0.D0
102   CONTINUE
      DO 103 IRAD=1,NRAD
        AREA(IRAD)=0.D0
103   CONTINUE
C
      IF (LEVGEO.EQ.1) THEN
C
C 1D SLAB-MODEL, DY = YDF, DZ = ZDF
C
        DO 110 IR=1,NR1STM
          AREA1(IR)=(RSURF(IR+1)-RSURF(IR))*YDF
          SX=(RSURF(IR+1)+RSURF(IR))*0.5
          IF (NLTRZ) THEN
            VOL(IR)=AREA1(IR)*ZDF
          ELSEIF (NLTRA) THEN
            VOL(IR)=AREA1(IR)*(SX+RMTOR)*TANAL/ALPHA*PI2A
          ELSE
            WRITE (iunout,*) 
     .        'INVALID OPTION IN SUBR. VOLUME, EXIT CALLED'
            CALL EXIT_OWN(1)
          ENDIF
110     CONTINUE
        GOTO 190
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
C 1D GRID OF CIRCLES, ELLIPSES OR TRIANGULAR FLUXSURFACES
C
        XNULL=0.
        IFLAG=1
        DNULL=0.
        DONE=1.
        CALL ARELLP(EP1(1),DNULL,ELL(1),DONE,TRI(1),DNULL,
     .              RSURF(1),DNULL,PI2A,XNULL,IFLAG,
     .              AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
        AREA1(0)=AELL
        DO 121 IR=1,NR1STM
C  AREA, CENTER OF GRAVITY
          IRP=IR+1
          CALL ARELLP(EP1(IRP),EP1(IR),ELL(IRP),ELL(IR),
     .                TRI(IRP),TRI(IR),
     .                RSURF(IRP),RSURF(IR),PI2A,XNULL,IFLAG,
     .                AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
          AREA1(IR)=AELL
          IF (NLTRT) THEN
            VOL(IR)=AREA1(IR)*(SX+RMTOR)*PI2A
          ELSEIF (NLTRZ) THEN
            VOL(IR)=AREA1(IR)*ZDF
          ELSEIF (NLTRA) THEN
            VOL(IR)=AREA1(IR)*(SX+RMTOR)*TANAL/ALPHA*PI2A
          ENDIF
121     CONTINUE
        GOTO 190
C
      ELSEIF (LEVGEO.EQ.3) THEN
C
C   1D GRID OF POLYGONS
C
        DO 139 K=1,NPPLG
          DO 139 J=NPOINT(1,K),NPOINT(2,K)-1
          AR=ARTRIA(0._DP,0._DP,XPOL(1,J),YPOL(1,J),
     .                          XPOL(1,J+1),YPOL(1,J+1))
          AREA1(0)=AREA1(0)+AR
139     CONTINUE
        AREA1(0)=ABS(AREA1(0))
C
        CALL ARPOLY(XPOL,YPOL,NRPLG,N1STS,1,NR1ST,AREAP,XCOM,YCOM)
        DO IR=1,NR1STM
          DO IP=1,NRPLG-1
            IN = IR + (IP-1)*NR1ST
            AREA(IN) = AREAP(IR,IP)
          END DO
        END DO
C
        DO 131 IR=1,NR1STM
          DO 132 K=1,NPPLG
            DO 132 J=NPOINT(1,K),NPOINT(2,K)-1
              AREA1(IR)=AREA1(IR)+AREAP(IR,J)
132       CONTINUE
131     CONTINUE
C
        IF (NLTRT) THEN
C   PARTICLES SEE A TORUS
          DO 134 IR=1,NR1STM
            XC=0.
            AR=0.
            DO 135 JP=1,NPPLG
              DO 135 J=NPOINT(1,JP),NPOINT(2,JP)-1
                IN = IR + (J-1)*NR1ST
                XC=XC+AREAP(IR,J)*XCOM(IN)
                AR=AR+AREAP(IR,J)
135         CONTINUE
            XC=XC/(AR+EPS60)
            VOL(IR)=AREA1(IR)*(XC+RMTOR)*PI2A
134       CONTINUE
        ELSEIF (NLTRZ) THEN
C   PARTICLES SEE A CYLINDER OF LENGTH DZ = ZDF
          DO 133 IR=1,NR1STM
            VOL(IR)=AREA1(IR)*ZDF
133       CONTINUE
C   PARTICLES SEE A TORUS APPROXIMATED BY NTTRAM STRAIGHT CYLINDERS
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 136 IR=1,NR1STM
            XC=0.
            AR=0.
            DO 137 JP=1,NPPLG
              DO 137 J=NPOINT(1,JP),NPOINT(2,JP)-1
                IN = IR + (J-1)*NR1ST
                XC=XC+AREAP(IR,J)*XCOM(IN)
                AR=AR+AREAP(IR,J)
137         CONTINUE
            XC=XC/(AR+EPS60)
            VOL(IR)=AREA1(IR)*(XC+RMTOR)*PI2AT
136       CONTINUE
        ENDIF
        GOTO 190
C
      ELSEIF (LEVGEO.EQ.4) THEN
C
C  GRID DEFINED BY FINITE ELEMENTS
C
        DO 150 IR=1,NTRII
          XCOM(IR) = (XTRIAN(NECKE(1,IR))+XTRIAN(NECKE(2,IR))+
     .                  XTRIAN(NECKE(3,IR)))/3.
          YCOM(IR) = (YTRIAN(NECKE(1,IR))+YTRIAN(NECKE(2,IR))+
     .                  YTRIAN(NECKE(3,IR)))/3.
          AR=0.5*(XTRIAN(NECKE(2,IR))*(YTRIAN(NECKE(3,IR))
     >           -YTRIAN(NECKE(1,IR)))+XTRIAN(NECKE(3,IR))*
     >           (YTRIAN(NECKE(1,IR))
     >           -YTRIAN(NECKE(2,IR)))+XTRIAN(NECKE(1,IR))*
     >           (YTRIAN(NECKE(2,IR))-YTRIAN(NECKE(3,IR))))
          AREA(IR) = AR
150     CONTINUE
C
C
C   PARTICLES SEE A TORUS
        IF (NLTRT) THEN
          DO 154 IR=1,NTRII
            VOL(IR)=AREA(IR)*(XCOM(IR)+RMTOR)*PI2A
154       CONTINUE
C   PARTICLES SEE A CYLINDER OF LENGTH DZ = ZDF
        ELSEIF (NLTRZ) THEN
          DO 153 IR=1,NTRII
            VOL(IR)=AREA(IR)*ZDF
153       CONTINUE
C   PARTICLES SEE A TORUS APPROXIMATED BY NTTRAM STRAIGHT CYLINDERS
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 155 IR=1,NTRII
            VOL(IR)=AREA(IR)*(XCOM(IR)+RMTOR)*PI2AT
155       CONTINUE
        ENDIF
C
      ELSEIF (LEVGEO.EQ.5) THEN
        TWOTHIRD=2.0D0/3.0D0
        DO ITET=1,NTET
C  CENTER OF MASS
          XTCEN(ITET)=0.D0
          YTCEN(ITET)=0.D0
          ZTCEN(ITET)=0.D0
          DO J=1,4
            IC=NTECK(J,ITET)
            XTCEN(ITET)=XTCEN(ITET) + XTETRA(IC)
            YTCEN(ITET)=YTCEN(ITET) + YTETRA(IC)
            ZTCEN(ITET)=ZTCEN(ITET) + ZTETRA(IC)
          END DO
          XTCEN(ITET)=XTCEN(ITET)*0.25D0
          YTCEN(ITET)=YTCEN(ITET)*0.25D0
          ZTCEN(ITET)=ZTCEN(ITET)*0.25D0

          IC1 = NTECK(1,ITET)
          IC2 = NTECK(2,ITET)
          IC3 = NTECK(3,ITET)
          IC4 = NTECK(4,ITET)
C  CALCULATE VOLUMES
          PC1(1:3)= (/ XTETRA(IC1), YTETRA(IC1), ZTETRA(IC1) /)
          PC2(1:3)= (/ XTETRA(IC2), YTETRA(IC2), ZTETRA(IC2) /)
          PC3(1:3)= (/ XTETRA(IC3), YTETRA(IC3), ZTETRA(IC3) /)
          PC4(1:3)= (/ XTETRA(IC4), YTETRA(IC4), ZTETRA(IC4) /)
          VOL(ITET) = CAL_VOL(PC1,PC2,PC3,PC4)
          IF (VOL(ITET) < -EPS10) THEN
            WRITE (iunout,*) ' WARNING ! '
            WRITE (iunout,*) ' VOL(',ITET,') < 0 VOL= ',VOL(ITET)
          END IF
          IF (SUM(NTBAR(1:4,ITET)) < 0) VOL(ITET) = 0._DP ! COLLAPSED TET
          VOL(ITET) = MAX(VOL(ITET),0._DP)
          AREA(ITET)=VOL(ITET)**TWOTHIRD
        END DO
C
      ELSEIF (LEVGEO.EQ.6) THEN
C
C  GENERAL GEOMETRY OPTION: PROVIDE CELL VOLUMES (CM**3)
C                      ON ARRAY VOL(IC),IC=1,NSURFM
C                     (ALSO PROVIDE CENTER OF CELL:
C                      XCOM(IC),YCOM(IC))
C
        CALL VOLUSR(NR1ST,VOL)
C
      ENDIF
C
190   CONTINUE
C
C
C
C  SET RADIAL SURFACE LABELING MESHES RHOSRF AND RHOZNE
C
      VOLSR=0.
      AREAR=0.
      IF (LEVGEO.EQ.1) THEN
C  RHOSRF(J) = RSURF(J)
        DO 195 J=1,NR1ST
          RHOSRF(J)=RSURF(J)
195     CONTINUE
      ELSEIF (LEVGEO.EQ.2) THEN
C  PI*RHOSRF(J)**2 = AREA ENCLOSED BY ORIGIN AND SURFACE J
        DO 196 J=1,NR1ST
          AREAR=AREAR+AREA1(J-1)
          RHOSRF(J)=SQRT(AREAR/PIA)
196     CONTINUE
      ELSEIF (LEVGEO.EQ.3) THEN
C  PI*RHOSRF(J)**2 = AREA ENCLOSED BY ORIGIN AND SURFACE J
        DO 197 J=1,NR1ST
          AREAR=AREAR+AREA1(J-1)
          RHOSRF(J)=SQRT(AREAR/PIA)
197     CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
C  RHOSRF AND RHOZNE ARE NOT DEFINED FOR FEM-OPTION
        RETURN
      ELSEIF (LEVGEO.EQ.5) THEN
C  RHOSRF AND RHOZNE ARE NOT DEFINED FOR TETRAHEDRON-OPTION
        RETURN
      ELSEIF (LEVGEO.EQ.6) THEN
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
        RETURN
      ENDIF
C
      DO 199 J=1,NR1STM
        RHOZNE(J)=0.5*(RHOSRF(J)+RHOSRF(J+1))
199   CONTINUE
C  MIRROR POINT
      RHOZNE(NR1ST)=RHOSRF(NR1ST)+0.5*(RHOSRF(NR1ST)-RHOSRF(NR1STM))
C
      IF (LEVGEO.LE.3) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'FLUX-SURFACE LABELING GRIDS'
        CALL MASRR2('  N, RHOSRF,RHOZNE    ',
     .                    RHOSRF,RHOZNE,NR1ST)
        CALL LEER(2)
      ENDIF
C
      RETURN
C
C  2D (R-THETA OR X-Y) VOLUME ELEMENTS
C
200   CONTINUE
C
      IF (LEVGEO.EQ.1) THEN
C
        KP=1
        DO 220 I=1,NR1STM
          NCELL1 = I+(      (KP-1)*NP2T3)*NR1P2
          VSAVE = VOL(NCELL1)
          DO 220 J=1,NP2NDM
            FAC2=(PSURF(J+1)-PSURF(J))/YDF
            NCELLJ=I+((J-1)+(KP-1)*NP2T3)*NR1P2
            XCOM(NCELLJ)=RHOZNE(I)
            YCOM(NCELLJ)=PHZONE(J)
            VOL(NCELLJ)=VSAVE*FAC2
220     CONTINUE
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
        IFLAG=1
        IT=1
        DO 240 IR=1,NR1STM
          IRP=IR+1
          DO 250 IP=1,NP2NDM
            IPP=IP+1
            CALL ARELLP(EP1(IRP),EP1(IR),ELL(IRP),ELL(IR),
     .                  TRI(IRP),TRI(IR),
     .                  RSURF(IRP),RSURF(IR),PSURF(IPP),PSURF(IP),IFLAG,
     .                  AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
C
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            AREA(NCELL)=AELL
            XCOM(NCELL)=SX
            YCOM(NCELL)=SY
250       CONTINUE
240     CONTINUE
C
        IF (NLTRT) THEN
          DO 262 I=1,NR1STM
            DO 262 J=1,NP2NDM
              K=1
              NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
              VOL(NCELL)=AREA(NCELL)*(XCOM(NCELL)+RMTOR)*PI2A
              IF (VOL(NCELL).GE.0.D0) GOTO 262
              WRITE (iunout,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
              CALL MASJ2('J,I             ',I,J)
C             CALL EXIT_OWN(1)
262       CONTINUE
        ELSEIF (NLTRZ) THEN
          DO 260 I=1,NR1STM
            DO 260 J=1,NP2NDM
              K=1
              NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
              VOL(NCELL)=AREA(NCELL)*ZDF
              IF (VOL(NCELL).GE.0.D0) GOTO 260
              WRITE (iunout,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
              CALL MASJ2('J,I             ',I,J)
C             CALL EXIT_OWN(1)
260       CONTINUE
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 261 I=1,NR1STM
            DO 261 J=1,NP2NDM
              K=1
              NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
              VOL(NCELL)=AREA(NCELL)*(XCOM(NCELL)+RMTOR)*PI2AT
              IF (VOL(NCELL).GE.0.D0) GOTO 261
              WRITE (iunout,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
              CALL MASJ2('J,I             ',I,J)
C             CALL EXIT_OWN(1)
261       CONTINUE
        ENDIF
C
      ELSEIF (LEVGEO.EQ.3) THEN
C
        IF (NLTRT) THEN
          DO 265 I=1,NR1STM
            DO 265 JP=1,NPPLG
              DO 265 J=NPOINT(1,JP),NPOINT(2,JP)-1
                K=1
                NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
                VOL(NCELL)=ABS(AREAP(I,J))*(XCOM(NCELL)+RMTOR)*PI2A
                IF (VOL(NCELL).GE.0.D0) GOTO 265
                WRITE (iunout,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
                CALL MASJ2('J,I             ',I,J)
C               CALL EXIT_OWN(1)
265       CONTINUE
        ELSEIF (NLTRA) THEN
          PI2AT=TANAL/ALPHA*PI2A
          DO 267 I=1,NR1STM
            DO 267 JP=1,NPPLG
              DO 267 J=NPOINT(1,JP),NPOINT(2,JP)-1
                K=1
                NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
                VOL(NCELL)=ABS(AREAP(I,J))*(XCOM(NCELL)+RMTOR)*PI2AT
                IF (VOL(NCELL).GE.0.D0) GOTO 267
                WRITE (iunout,*) 'ERROR IN SUBR. VOLUME, VOL.LT.0'
                CALL MASJ2('J,I             ',I,J)
C               CALL EXIT_OWN(1)
267       CONTINUE
        ELSEIF (NLTRZ) THEN
          DO 268 I=1,NR1STM
            DO 268 JP=1,NPPLG
              DO 268 J=NPOINT(1,JP),NPOINT(2,JP)-1
                K=1
                NCELL=I+((J-1)+(K-1)*NP2T3)*NR1P2
                VOL(NCELL)=ABS(AREAP(I,J))*ZDF
268       CONTINUE
        ENDIF
C
      ELSEIF (LEVGEO.EQ.4) THEN
C
C  FINITE ELEMENT OPTION: NOTHING TO BE DONE HERE
C
      ELSEIF (LEVGEO.EQ.5) THEN
C
C  TETRAHEDRON OPTION: NOTHING TO BE DONE HERE
C
      ELSEIF (LEVGEO.EQ.6) THEN
C
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
C
      ENDIF
C
      RETURN
C
C    2D (R-Z), (R-PHI) OR (X-Z) VOLUME ELEMENTS
C
300   CONTINUE
C
      IF (NLTRZ.OR.NLTRA.OR.NLTRT) THEN
C
        IT=1
        DO 320 J=1,NP2ND
        DO 320 I=1,NR1ST
          NCELL1 = I+((J-1)            )*NR1P2
          VSAVE=VOL(NCELL1)
          DO 320 K=1,NT3RDM
            FAC3=(ZSURF(K+1)-ZSURF(K))/ZDF
            NCELLK=I+((J-1)+(K-1)*NP2T3)*NR1P2
            VOL(NCELLK)=VSAVE*FAC3
320     CONTINUE
C
      ENDIF
C
      RETURN
C
C   ADDITIONAL CELL VOLUMES
C
400   CONTINUE
C
C   ADDITIONAL CELL VOLUMES ARE DEFAULTED TO 1. AT PRESENT
C
      DO 410 J=NSURF+1,NSBOX
        VOL(J)=1.
410   CONTINUE
      IF (NLADD) THEN
        DO 411 J=NSURF+1,NSBOX
          IF (VOLADD(J-NSURF).GT.0.D0) VOL(J)=VOLADD(J-NSURF)
411     CONTINUE
      ENDIF

      DEALLOCATE (AREAP)
c slmod begin - debug
c      DEALLOCATE (AREA1)
c slmod end
C
      RETURN
C
999   CONTINUE
      WRITE (iunout,*) 'UNWRITTEN OPTION CALLED IN SUBR. VOLUME '
      CALL EXIT_OWN(1)
      END
