C
!pb  11.12.06:  allow letters 'f' or 't' in case name of fem or tetrahedron 
!pb             calculation
!pb  27.12.06:  bug fix: increase NSTS in case of time dependent mode
!pb  15.01.07:  additional line in input block 4 defining HYDKIN default model
!pb  02.03.07:  NUMSEC=4 introduced
!pb  20.03.07:  include input block written by HYDKIN default model
!pb  22.03.07:  input for NLFEM and NLTET corrected. 
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
      INTEGER, ALLOCATABLE :: INDSRC(:), IEIGEN(:)
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
     .           NS2, NS3, INM1, INM2, INM3, INMDL, IEND, ITOK, IER,
     .           N_REAC, N_SPEC, N_ATOMS, N_MOL, N_IONS, N_TESTIONS,
     .           N_BULKIONS, NB4, NS4, INM4, IUNIN_SAVE, I1
      REAL(DP) :: SORIND, SORLIM, DUMM1, ROA, ZAA, ZZA, ZGA, YAA, YYA,
     .            ZIA, YP, XP, YIA, YGA
      LOGICAL :: NLSCL, NLTEST, NLANA, NLDRFT, NLCRR, NLERG, NLIDENT,
     .           LHABER, NLONE, LTSTV, LINCLUDE
      LOGICAL :: NLSLB, NLCRC,  NLELL, NLTRI,  NLPLG, NLFEM, NLTET,
     .           NLGEN
      LOGICAL :: NLRAD, NLPOL,  NLTOR, NLADD,  NLMLT, NLTRIM
      LOGICAL :: NLTRA, NLTRT, NLTRZ
      LOGICAL :: PLTL2D, PLTL3D, LRPSCUT, LHYDDEF
      CHARACTER(80) :: ZEILE, CASENAME, FILENAME, ULINE, FILE
      CHARACTER(12) :: HYDKIN_DEFAULT, CHR
      CHARACTER(4) :: CLAB
      CHARACTER(15), ALLOCATABLE :: HYDSPEC(:)
      CHARACTER(1000) :: HLINE
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
      NREAC_LINES=0
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
     .                    NLERG,NLIDENT,NLONE,LTSTV
      ELSE
        READ (ZEILE,6665) NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
     .                    NLERG,NLIDENT,NLONE,LTSTV
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
        READ (IUNIN,6666) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
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

          IF (NLFEM .OR. NLTET) THEN  
            READ (IUNIN,'(A72)') ZEILE
            CLAB = ZEILE(1:4)
            CALL UPPERCASE(CLAB)
            IF (INDEX(ZEILE,'CASE')==0) THEN
              WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION '
              WRITE (IUNOUT,*) 
     .          ' TRIANGLE OR TETRAHEDRON GRID SWITCHED ON '
              WRITE (IUNOUT,*) ' BUT NO CASENAME SPECIFIED '
              CALL EXIT_OWN(1)
            END IF

            READ (ZEILE(6:),'(A66)') CASENAME
            CASENAME=ADJUSTL(CASENAME)
            LL=LEN_TRIM(CASENAME)
          END IF  
          
          IF (NLFEM) THEN
            NTRII=NR1ST
            NTRI=MAX(NTRI,NR1ST)
            
            FILENAME=CASENAME(1:LL) // '.npco_char'
            OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')
            ZEILE='*   '      
            DO WHILE (ZEILE(1:1) == '*')
              READ (30,'(A)') ZEILE
            END DO

            READ (ZEILE,*) NRKNOT
            CLOSE (UNIT=30)
            NKNOT=NRKNOT

            FILENAME=CASENAME(1:LL) // '.elemente'
            OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')

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
     .            FORM='FORMATTED')
              
            ZEILE='*   '
            DO WHILE (ZEILE(1:1) == '*')
               READ (30,'(A100)') ZEILE
            END DO
              
            NGITT=1
            DO I=1,NTRII
              READ (30,*) ID, NB1, NS1, INM1,
     .                        NB2, NS2, INM2, 
     .                        NB3, NS3, INM3
              IF (INM1 /= 0) NGITT = NGITT + 1
              IF (INM2 /= 0) NGITT = NGITT + 1
              IF (INM3 /= 0) NGITT = NGITT + 1
            END DO
              
            CLOSE (UNIT=30)
              
          ENDIF
          
          IF (NLTET) THEN
            NTET=NR1ST
            NTETRA=MAX(NTETRA,NR1ST)
            
            FILENAME=CASENAME(1:LL) // '.npco_char'
            OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')
            ZEILE='*   '      
            DO WHILE (ZEILE(1:1) == '*')
              READ (30,'(A)') ZEILE
            END DO

            READ (ZEILE,*) NCOOR
            CLOSE (UNIT=30)
            NCOORD=NCOOR

            FILENAME=CASENAME(1:LL) // '.elemente'
            OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')

            ZEILE='*   '      
            DO WHILE (ZEILE(1:1) == '*')
              READ (30,'(A)') ZEILE
            END DO

            READ (ZEILE,*) NTET
            CLOSE (UNIT=30)
            NTETRA=NTET+1
            NR1ST=NTETRA

            FILENAME=CASENAME(1:LL) // '.neighbors'
            OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .            FORM='FORMATTED')
              
            ZEILE='*   '
            DO WHILE (ZEILE(1:1) == '*')
               READ (30,'(A100)') ZEILE
            END DO
              
            NGITT=1
            DO I=1,NTET
              READ (30,*) ID, NB1, NS1, INM1,
     .                        NB2, NS2, INM2, 
     .                        NB3, NS3, INM3,
     .                        NB4, NS4, INM4
              IF (INM1 /= 0) NGITT = NGITT + 1
              IF (INM2 /= 0) NGITT = NGITT + 1
              IF (INM3 /= 0) NGITT = NGITT + 1
              IF (INM4 /= 0) NGITT = NGITT + 1
            END DO
              
            CLOSE (UNIT=30)
              
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

! CHECK FOR INCLUDE LINE
      READ (IUNIN,'(A72)') ZEILE
      ULINE = ZEILE
      IREAD=1
      CALL UPPERCASE(ULINE)
      I1 = INDEX(ULINE,'INCLUDE')
      LINCLUDE = .FALSE.
      IF (I1 > 0) THEN
        IREAD = 0
        CALL READ_TOKEN(ZEILE(I1+7:),' ',FILE,ITOK,IER,.FALSE.)
        LINCLUDE = .TRUE.
        IUNIN_SAVE = IUNIN
        IUNIN = 2
        OPEN (IUNIN,FILE=FILE,FORM='FORMATTED',ACCESS='SEQUENTIAL')
      END IF          
C
C
      IF (IREAD == 0) READ (IUNIN,*)
      WRITE (iunout,*) 
     .  '       ATOMIC REACTION CARDS, NREACI DATA FIELDS'
      READ (IUNIN,'(A72)') ZEILE
      CALL UPPERCASE(ZEILE)
      IEND=INDEX(ZEILE,'DEFAULT') 
      LHYDDEF =.FALSE.
      IF (IEND > 0) THEN
        CALL READ_TOKEN(ZEILE(IEND+7:),' ',HYDKIN_DEFAULT,ITOK,IER,
     .                  .FALSE.)
        LHYDDEF=.TRUE.
        READ (IUNIN,'(A72)') ZEILE
      END IF
      READ (ZEILE,*) NREACI
      NREAC = MAX(NREAC,NREACI)+NREAC_ADD
C
      NREAC_LINES=0
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        NREAC_LINES=NREAC_LINES+1
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
     .        INDEX(ULINE,'BOLTZMANN')+INDEX(ULINE,'COLRAD')+
     .        INDEX(ULINE,'CONSTANT')
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

      IF (LHYDDEF) THEN
        HYDKIN_DEFAULT=ADJUSTL(HYDKIN_DEFAULT)
        LL=LEN_TRIM(HYDKIN_DEFAULT)
        FILENAME=HYDKIN_DEFAULT(1:LL) // '.reactions'
        OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED')
        READ (30,*)
        READ (30,*) CHR,n_reac
        READ (30,*) CHR,n_spec
        READ (30,*) CHR,n_atoms
        READ (30,*) CHR,n_ions
        READ (30,*) CHR,n_mol

        ALLOCATE (HYDSPEC(N_SPEC))
        ALLOCATE (IEIGEN(N_SPEC))

        READ (30,*) HYDSPEC(1:N_SPEC)
        READ (30,'(A1000)') HLINE
        READ (HLINE(52:),*) IEIGEN(1:N_SPEC)

!pb        N_BULKIONS = COUNT((SCAN(HYDSPEC(1:N_SPEC),'+') > 0) .AND. 
!pb     .                     (IEIGEN(1:N_SPEC) == 0))
!pb        N_TESTIONS = N_IONS - N_BULKIONS

        N_BULKIONS = COUNT(IEIGEN(1:N_SPEC) == 0)
        N_TESTIONS = COUNT((SCAN(HYDSPEC(1:N_SPEC),'+-') > 0) .AND. 
     .                     (IEIGEN(1:N_SPEC) /= 0))

        NATM = NATM + N_ATOMS
        NMOL = NMOL + N_MOL
        NION = NION + N_TESTIONS
        NPLS = NPLS + N_BULKIONS
        NREAC = NREAC + N_REAC
        NREAC_LINES = NREAC_LINES + N_REAC
      END IF
        

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

      IF (LINCLUDE) THEN
        CLOSE (IUNIN)
        IUNIN = IUNIN_SAVE
        LINCLUDE =.FALSE.
      END IF

      DO 
        READ (IUNIN,'(A72)') ZEILE
        IF ((ZEILE(1:3) == '***') .AND. 
     ,      (INDEX(ZEILE,'6.') > 0)) EXIT
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
          IF (NTIME >= 1) NSTS=NSTS+1
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

      END
