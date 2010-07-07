cc  coupling to fidap: ein brick = 48 tetraeder, so dass
cc  mindestens jeder fidap knoten eine ecke eines tetraeder ist.
cc  durch die zerlegung in 48 (nicht nur 5 oder 6) tetraeder koennen
cc  die bricks (quader) beliebig orientiert sein, es geht immer so auf,
cc  dass eine brick-seitenflaeche, zerlegt in tetraeder-dreieckseiten,
cc  auf die zerlegung der nachbarflaeche passt.
cc
C   THIS CODE SEGMENT CONTAINES VARIOUS SUBROUTINES NEEDED FOR
C   INTERFACING THE EIRENE CODE TO PLASMA FLUID CODES.
C   IT READS GEOMETRICAL DATA (MESHES) FROM FILE FT30
C   AND PRODUCES THE EIRENE INPUT DATA (BLOCK 2).
C   IT READS PLASMA BACKGROUND DATA FROM FILE (FT31) OR COMMON BLOCKS,
C   IT THEN PRODUCES INPUT DATA FOR EIRENE
C   INPUT BLOCK 5 (PLASMA DATA) AND BLOCK 7 (SURFACE RECYCLING SOURCES)
C
C
*DK COUPLE
C
      SUBROUTINE INFCOP
C
C
C     THIS SUBROUTINE DEFINES THE PLASMA MODEL IN CASE OF A COUPLED
C     NEUTRAL-PLASMA CALCULATION
C
C     THE ENTRY "IF0COP" RECEIVES GEOMETRICAL INPUT DATA FROM AN
C     EXTERNAL FILE (E.G. OTHER PLASMA CODES)
C     AND PREPARES THEM FOR AN EIRENE RUN
C
C     THE ENTRY "IF1COP" RECEIVES PLASMA INPUT DATA FROM AN
C     EXTERNAL FILE (E.G. OTHER PLASMA CODES)
C     AND PREPARES THEM FOR AN EIRENE RUN
C
C     THE ENTRY "IF2COP" PREPARES THE SOURCE SAMPLING DISTRIBUTION
C     FROM THE EXTERNAL DATA, AND MAY OVERWRITE OTHER INPUT
C     DATA FROM BLOCKS 1 TO 13 AS WELL
C
C     THE ENTRIES "IF3COP, IF4COP" RETURN  RESULTS TO AN EXTERNAL FILE
C
C
C
      USE PRECISION
      USE PARMMOD
      USE BRASPOI
      USE COMUSR
      USE CESTIM
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CPLOT
      USE CINIT
      USE CPOLYG
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CCOUPL
      USE CGEOM
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE CTETRA
      USE COMPRT
      USE COMNNL
      USE COMSOU
      USE CSTEP
      USE CTEXT
      USE CLGIN
      USE COUTAU
      USE COMXS
      USE CSPEI
      USE CTRIG
      USE module_avltree

      IMPLICIT NONE

      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL
C
      REAL(DP) :: SCALN(0:NFL)
      REAL(DP) :: DI(NPLS), VP(NPLS), CSP(NPLS)
      REAL(DP) :: EPEL(NRTAL), HELPAR(NRTAL), DIVI(NRTAL)
      REAL(DP) :: EMSAMP(NCOOR), EM1SMP(NCOOR), EM2SMP(NCOOR)
      REAL(DP) :: ABSAMP(NCOOR), AB1SMP(NCOOR), AB2SMP(NCOOR)
      REAL(DP) :: RIX(NCOOR), RIY(NCOOR), RIZ(NCOOR)
      REAL(DP) :: XKAP(NCOOR), YKAP(NCOOR), ZKAP(NCOOR)
      REAL(DP) :: EMIS(NCOOR), ABSORB(NCOOR), VOLSUM(NCOOR)
C
      REAL(DP) :: PUX(NRAD),PUY(NRAD),PVX(NRAD),PVY(NRAD),
     R          TORL(NSTRA,NGITT),EFLX(NSTRA),
     R          DUMMY(0:NDXP,0:NDYP),
     R          ESHT(NSTEP,NGITT),ORI(NSTEP,NGITT),
     R          SFNIT(0:NSTEP),SFEIT(0:NSTEP),
     R          SFEET(0:NSTEP),SHEAE(0:NSTEP),SHEAI(0:NSTEP)
      INTEGER, ALLOCATABLE :: NRWL(:)
C
      LOGICAL :: LSHORT=.FALSE., LSTOP=.TRUE., LSTP, LSHIFT

      INTEGER :: NBRICK, NQUAD, NLINE
      COMMON /CINTFA/ NBRICK, NQUAD, NLINE

      INTEGER, ALLOCATABLE, SAVE :: NODBRICK(:,:), NODQUAD(:,:)

      INTEGER, ALLOCATABLE :: IHELP(:), NOSTS(:), NQDSTS(:), NOSPEC(:),
     .                        NOSTEP(:), NQDSTP(:)
      REAL(DP) :: TIS(NPLS),DIS(NPLS),VXS(NPLS),VYS(NPLS),VZS(NPLS)
      REAL(SP) :: PARM(16)
      INTEGER :: ITSIDE(3,4)

      REAL(DP) :: EMAXW, ESUM, ETOT, EADD, PERWI, PARWI, FLX, FLXI, 
     .          VPX, VPY, GAMMA, SHEATH, VT, PN1, VTEST, VR, DRR,
     .          PERW, PARW, CS, CUR, X1, Y1, Z1, X2, Y2, Z2, PHI1, PHI2,
     .          STEP, THMAX, BXS, BYS, BZS, PM1, TE, VPZ, ZMFPTH,
     .          EEMAX, RP1, OR, EESHT, V, T, TES, VL, XCO, YCO, ZCO,
     .          DIST, DIST1, DIST2, DIST3, DIST4, BVAC, RECADD, EEADD,
     .          SUMN, SUMEE, RECTOT, FTABRC1, FEELRC1, XGAU, YGAU, ZGAU,
     .          DIXDX, DIYDY, DIZDZ, DETI, FAC, SIGMX
      REAL(DP) :: TQVEC(27), XVEC(27), YVEC(27), ZVEC(27), GK(3), 
     .          DFVEC(3,27), AJ(3,3), AJM1(3,3), RIXVEC(27), RIYVEC(27), 
     .          RIZVEC(27), E(3,3), AJM1T(3,3), TEVEC(27), 
     .          DF2VEC(3,3,27), chrsym(3,3,3)
      real(dp) :: dtdxyz(3,27)
      real(dp) :: u1(3), u2(3), u3(3), eu1(3),eu2(3),eu3(3),guik(3,3),
     .            gradphi(3), eo1(3), eo2(3), eo3(3), gradphi2(3), 
     .            goik(3,3), dfguik(3,3,3)
      real(dp) :: dxdr, dxds, dxdt, dydr, dyds, dydt, dzdr, dzds, dzdt,
     .            e123, dphidx, dphidy, dphidz, weifakt, sumik,
     .            xlaplace, test, testte, dtedxj,flaplace
      REAL(DP) :: TEG, DTEDXG, DTEDYG, DTEDZG, gikd2te, sumj, suml
      INTEGER :: IFIRST, NREC11, NDXY, NR1STQ, ISTRAI, IY, IRC, NFALSE,
     .           JC, J, NDXYM, NDYAM, NDXAM, N1, N2, N3IN, N3EN, ITFAL,
     .           IAOT, IAIN, IR, IP, IR1, IP1, K, MPER, N3, NT, IERROR,
     .           IPL, IMODE, IFRST, ITARG, INEW, IOLD, ISP, ISTEP,
     .           N1EN, N1IN, N2EN, N2IN, IPRT, NDZA, I, NTGPRI, IT, 
     .           IEPLS, I4, IIPLS, IG, NPES, ITRI, NTACT, IN, NTOLD, IS,
     .           ISTS, IC, IDU, I1, I2, I3, NOD, NTT, IEL, NBR, M,
     .           NQU, IQ, IMATCH, LANF, LEND, ICO, ITET, ICOLUMN, IANF,
     .           NREAD, ISRFSI, ISR, IIRC, IRRC, INC, IPO, L, II, JJ,
     .           IPLSTI, IPLSV, IPLV
      INTEGER :: INDCO(27),INDF1(9),INDF2(9),INDF3(9),INDF4(9),INDF5(9),
     .           INDF6(9),NSORQUAD(9)
      INTEGER, INTENT(IN) :: ISTRAA, ISTRAE, NEW_ITER
      INTEGER, EXTERNAL :: IDEZ
      REAL(DP), ALLOCATABLE :: TEF(:),DENF(:,:)
      REAL(DP), ALLOCATABLE, SAVE :: DRSTDXYZ(:,:,:),
     .                         DTEDX(:),DTEDY(:),DTEDZ(:),XLAPLA(:)
C
      CHARACTER(72) :: ZEILE, NAMENT, SPECNAME
      CHARACTER(10) :: CHR, FORM
      CHARACTER(6) :: CITARG
      CHARACTER(72), ALLOCATABLE :: ENTITY(:)
      CHARACTER(1000) :: LINE
      LOGICAL, ALLOCATABLE :: LUSED(:)
C
C
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/
C
C
      ENTRY IF0COP
C
      LSHORT=.FALSE.
C
      GOTO 99990
C
C  TO INITALIZE THE SHORT CYCLING, THE GEOMETRY HAS TO BE
C  DEFINED ONCE (ENTRY: INTER0)
C
      ENTRY INTER0
      LSHORT=.TRUE.
99990 CONTINUE
C
      IERROR=0
C
      IMODE=IABS(NMODE)
C
      IF (.NOT.LSHORT.AND.ITIMV.LE.1) THEN
        WRITE (iunout,*) '        SUBROUTINE INFCOP IS CALLED  '
C  READ INPUT DATA OF BLOCK 14
C  SAVE INPUT DATA OF BLOCK 14 FOR SHORT CYCLE ON COMMON CCOUPL
        CALL LEER(1)
        CALL ALLOC_CCOUPL(1)
        READ (IUNIN,'(5L1)') LSYMET,LBALAN
        IF (TRCINT)
     .  WRITE (iunout,*) ' LSYMET,LBALAN = ',LSYMET,LBALAN
        DO 20 IPL=1,NPLSI
          READ (IUNIN,'(2I6,2E12.4)') I,IFLB(IPL),FCTE(IPL),BMASS(IPL)
          IF (TRCINT)
     .    WRITE (iunout,*) IPL,IFLB(IPL),FCTE(IPL),BMASS(IPL)
20      CONTINUE
C  NUMBER OF DIFFERENT ENTITIES: NTARGI
        READ (IUNIN,'(I6)') NTARGI
        WRITE (iunout,*) '        NTARGI= ',NTARGI
        CALL LEER(1)
        ALLOCATE (NOSTS(NTARGI))
        ALLOCATE (NOSTEP(NTARGI))
        ALLOCATE (ENTITY(NTARGI))
        DO 30 IT=1,NTARGI
          nament=repeat(' ',72)
          READ (IUNIN,'(3I6,1X,A72)') I,NOSTS(IT),NOSTEP(IT),NAMENT
          LANF=verify(nament,' ')
          LEND=verify(nament,' ',.true.)
          ENTITY(IT)=repeat(' ',72)
          ENTITY(IT)(1:lend-lanf+1) = nament(lanf:lend)
          WRITE (iunout,*) ' ENTITY ',ENTITY(IT)
          CALL UPPERCASE (ENTITY(IT))
          WRITE (iunout,*) ' ENTITY ',ENTITY(IT)
          IF ((NOSTEP(IT) < 0).OR.(NOSTEP(IT) > NTARGI)) THEN
            WRITE (iunout,*) ' NOSTEP MUST BE >= 0 AND <= NTARGI '
            WRITE (iunout,*) ' NOSTEP, NTARGI ', NOSTEP(IT), NTARGI 
            CALL EXIT(1)
          END IF
          IF (TRCINT)
     .      WRITE (iunout,'(3I6,1X,A72)') 
     .            IT,NOSTS(IT),NOSTEP(IT),ENTITY(IT)
          IF (TRCINT) CALL LEER(1)
30      CONTINUE
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM FIDAP INTO EIRENE
C  HERE: FIDAP VOLUME TALLIES
        READ (IUNIN,'(I6)') NAINB
        NAIN = MAX(NAIN,NAINB)
        CALL ALLOC_CCOUPL(2)
        WRITE (iunout,*) '        NAINI = ',NAINB
        IF (NAINB.GT.NAIN) THEN
          CALL MASPRM ('NAIN',4,NAIN,'NAINB',5,NAINB,IERROR)
          WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
          CALL EXIT(1)
        ENDIF
        IF (TRCINT.AND.NAINB.GT.0)
     .      WRITE (iunout,*) 'I,NAINS(IAIN),NAINT(IAIN)'
        DO 40 IAIN=1,NAINB
          READ (IUNIN,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
          READ (IUNIN,'(A72)') TXTPLS(IAIN,12)
          READ (IUNIN,'(2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
          IF (TRCINT) THEN
            WRITE (iunout,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
            WRITE (iunout,'(1X,A72)') TXTPLS(IAIN,12)
            WRITE (iunout,'(1X,2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
          ENDIF
40      CONTINUE
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO FIDAP
C  HERE: EIRENE SURFACE TALLIES
        READ (IUNIN,'(I6)') NAOTB
        WRITE (iunout,*) '        NAOTI = ',NAOTB
        IF (NAOTB.GT.NLIMPS) THEN
          CALL MASPRM ('NLIMPS',6,NLIMPS,'NAOTB',5,NAOTB,IERROR)
          WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
          CALL EXIT(1)
        ENDIF
        IF (TRCINT.AND.NAOTB.GT.0)
     .      WRITE (iunout,*) 'I,NAOTS(IAOT),NAOTT(IAOT)'
        DO 50 IAOT=1,NAOTB
          READ (IUNIN,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
          IF (TRCINT) THEN
            WRITE (iunout,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
          ENDIF
50      CONTINUE

        READ (IUNIN,'(6E12.4)') ZMFPTHI, SIGMX, TDGTEMX
      ENDIF
C
C READING BLOCK 14 FROM FORMATTED INPUT FILE (IUNIN) FINISHED
C
C SAVE SOME MORE INPUT DATA FOR SHORT CYCLE ON COMMON CCOUPL
      LNLPLG=NLPLG
      LNLDRF=NLDRFT
      LTRCFL=TRCFLE
      NSTRI=NSTRAI
      DO 60 ISTRA=1,NSTRAI
        LNLVOL(ISTRA)=NLVOL(ISTRA)
60    CONTINUE
      NMODEI=NMODE
      NFILNN=NFILEN
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                              AND IN SUBR. COLLIDE)
      NCPVI=max(NPLSI,13)
      IF (NCPVI.GT.NCPV) THEN
        WRITE (iunout,*) 'FROM INTERFACING SUBROUTINE INFCOP: '
        CALL MASPRM('NCPV',4,NCPV,'NCPVI',5,NCPVI,IERROR)
        CALL EXIT(1)
      ENDIF

      ICPVE(1)=3
      ICPRC(1)=1
      TXTTAL(1,NTALM)=
     .  'TOTAL SAMPLED EMISSION (copv(1)=eppht ?)                    '
      TXTSPC(1,NTALM)='                          '
      TXTUNT(1,NTALM)='W/CM**3                   '

      ICPVE(2)=3
      ICPRC(2)=1
      TXTTAL(2,NTALM)=
     .  'TOTAL SAMPLED ABSORPTION (copv(2))                          '
      TXTSPC(2,NTALM)='                          '
      TXTUNT(2,NTALM)='W/CM**3                   '

      ICPVE(3)=3
      ICPRC(3)=1
      TXTTAL(3,NTALM)=
     .  'EMISSION, THICK FRACTION (copv(3))                          '
      TXTSPC(3,NTALM)='                          '
      TXTUNT(3,NTALM)='W/CM**3                   '

      ICPVE(4)=3
      ICPRC(4)=1
      TXTTAL(4,NTALM)=
     .  'ABSORPTION THICK FRACTION (copv(4))                         '
      TXTSPC(4,NTALM)='                          '
      TXTUNT(4,NTALM)='W/CM**3                   '

      ICPVE(5)=3
      ICPRC(5)=1
      TXTTAL(5,NTALM)=
     .  'EMISSION THIN FRACTION (copv(5))                            '
      TXTSPC(5,NTALM)='                          '
      TXTUNT(5,NTALM)='W/CM**3                   '

      ICPVE(6)=3
      ICPRC(6)=1
      TXTTAL(6,NTALM)=
     .  'ABSORPTION THIN FRACTION (copv(6))                          '
      TXTSPC(6,NTALM)='                          '
      TXTUNT(6,NTALM)='W/CM**3                   '

      ICPVE(7)=1
      ICPRC(7)=1
      TXTTAL(7,NTALM)=
     .  'I_x (copv(7))                                               '
      TXTSPC(7,NTALM)='                          '
      TXTUNT(7,NTALM)='AMP                       '

      ICPVE(8)=1
      ICPRC(8)=1
      TXTTAL(8,NTALM)=
     .  'I_y (copv(8))                                               '
      TXTSPC(8,NTALM)='                          '
      TXTUNT(8,NTALM)='AMP                       '

      ICPVE(9)=1
      ICPRC(9)=1
      TXTTAL(9,NTALM)=
     .  'I_z (copv(9))                                               '
      TXTSPC(9,NTALM)='                          '
      TXTUNT(9,NTALM)='AMP                       '

      ICPVE(10)=0
      ICPRC(10)=0
      TXTTAL(10,NTALM)=
     .  'Kappa   (copv(10))                                          '
      TXTSPC(10,NTALM)='                          '
      TXTUNT(10,NTALM)='AMP                       '

      ICPVE(11)=0
      ICPRC(11)=0
      TXTTAL(11,NTALM)=
     .  'Laplace Operator (copv(11))                                 '
      TXTSPC(11,NTALM)='                          '
      TXTUNT(11,NTALM)='AMP                       '

      ICPVE(12)=3
      ICPRC(12)=1
      TXTTAL(12,NTALM)=
     .  'TOTAL SAMPLED ABSORPTION BY COLLISION ESTIMATOR (copv(12))  '
      TXTSPC(12,NTALM)='                          '
      TXTUNT(12,NTALM)='W/CM**3                   '

      ICPVE(13)=3
      ICPRC(13)=1
      TXTTAL(13,NTALM)=
     .  'TOTAL SAMPLED ABSORPTION BY TRACKLENGTH ESTIMATOR (copv(13))'
      TXTSPC(13,NTALM)='                          '
      TXTUNT(13,NTALM)='W/CM**3                   '
C
      OPEN (UNIT=29,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 29
C
      OPEN (UNIT=30,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 30
C
C  READ IN DATA TO SET UP GEOMETRY FOR NEUTRAL GAS TRANSPORT CODE
C  STATEMENT NUMBER 1000 ---> 1999
C
C  AT PRESENT THE DATA COME FROM THE FILE FT30
C  THIS PART WILL HAVE TO BE MODIFIED AS SOON AS BRAAMS PROVIDES
C  CELL VERTICES AND CUT DESCRIBTION
C
1000  CONTINUE
C
      INMP1I=0
      INMP2I=0
      INMP3I=0
      NTET = 0
      NCOOR = 0
      NTBAR = 0
      NTSEITE = 0
      ncltet = 0
      itethand=1
      ALLOCATE (COORTET(NCOORD))
      DO I=1,NCOORD
        NULLIFY(COORTET(I)%PTET)
      END DO

      nr1p2 = nr1st
      np2t3 = np2nd
      nr1ori = nr1st
      np2ori = np2nd
      nt3ori = nt3rd
      IF (SUM(NOSTEP(1:NTARGI)) > 0) THEN
        NGITT = NQUAD*8
      ELSE
        NGITT = 1
      END IF

      CALL ALLOC_CSTEP

      RRSTEP(1:NSTEP,1) = 0.D0
      ALLOCATE(NRWL(NTARGI))
      NRWL=0

      READ (30,'(A72)') ZEILE
      DO WHILE (INDEX(ZEILE,'NODAL COORDINATES') == 0)
         READ (30,'(A72)') ZEILE
      END DO

      DO IC=1,NCOORD
         READ (30,*) IDU,XTETRA(IC),YTETRA(IC),ZTETRA(IC)
      END DO
      NCOOR=NCOORD

!  CONVERT TO CM 
      XTETRA(1:NCOORD)=XTETRA(1:NCOORD)*0.1_DP
      YTETRA(1:NCOORD)=YTETRA(1:NCOORD)*0.1_DP
      ZTETRA(1:NCOORD)=ZTETRA(1:NCOORD)*0.1_DP

      write (iunout,*) ' xmin = ',minval(xtetra(1:ncoor))
      write (iunout,*) ' xmax = ',maxval(xtetra(1:ncoor))
      write (iunout,*) ' ymin = ',minval(ytetra(1:ncoor))
      write (iunout,*) ' ymax = ',maxval(ytetra(1:ncoor))
      write (iunout,*) ' zmin = ',minval(ztetra(1:ncoor))
      write (iunout,*) ' zmax = ',maxval(ztetra(1:ncoor))

      ALLOCATE (NODBRICK(27,NBRICK))
      ALLOCATE (NODQUAD(9,NQUAD))
      ALLOCATE (NQDSTS(NQUAD))
      ALLOCATE (NQDSTP(NQUAD))
      ALLOCATE (LUSED(NQUAD))
      NODBRICK = 0
      NODQUAD = 0
      NQDSTS = 0
      NQDSTP = 0
      LUSED = .FALSE.

      NBR = 0
      NQU = 0

      READ (30,'(A72)') ZEILE
      DO WHILE (INDEX(ZEILE,'ELEMENT GROUPS') == 0)
         READ (30,'(A72)') ZEILE
      END DO

      DO
        READ (30,'(A72)',END=9) ZEILE
        I1 = INDEX(ZEILE,'ELEMENTS:')
        DO WHILE (I1 == 0)
          READ (30,'(A72)',END=9) ZEILE
          I1 = INDEX(ZEILE,'ELEMENTS:')
        END DO
        I2 = INDEX(ZEILE,'NODES:')
        I3 = INDEX(ZEILE,'GEOMETRY:')
        READ (ZEILE(I2+7:I3-1),*) NOD

        READ (ZEILE(I1+10:I2),*) NTT

        READ (30,'(A72)') ZEILE
        I1 = INDEX(ZEILE,'ENTITY NAME:')
        DO WHILE (I1 == 0)
          READ (30,'(A72)') ZEILE
          I1 = INDEX(ZEILE,'ENTITY NAME:')
        END DO
        
        LANF=i1+12+verify(zeile(i1+12:),' ')-1
        LEND=verify(zeile,' ',.true.)
        NAMENT = ZEILE(lanf:lend)
        write (iunout,*) ' nament: ',nament
        CALL UPPERCASE(NAMENT)
        ISTS=0
        ISTEP=0
        DO IT=1,NTARGI
          IF (ENTITY(IT) == NAMENT) THEN
            ISTS=NOSTS(IT)
            ISTEP=NOSTEP(IT)
          END IF
        END DO
        
        IF (NOD == 27) THEN
! READ COORDINATE NUMBERS OF BRICKS
          DO IEL=1,NTT
            NBR = NBR + 1
            READ (30,*) IDU,(NODBRICK(J,NBR),J=1,NOD)
          END DO
        ELSE IF (NOD == 9) THEN
! READ COORDINATE NUMBERS OF QUADRILATERALS
          DO IEL=1,NTT
            NQU = NQU + 1
            READ (30,*) IDU,(NODQUAD(J,NQU),J=1,NOD)
!pb            CALL BUBBLE (NODQUAD(1:NOD,NQU),NOD)
            NQDSTS(NQU) = NLIM+ISTS
            NQDSTP(NQU) = ISTEP
          END DO
        ELSE
          EXIT
        END IF
         
      END DO

 9    CONTINUE

      NTACT=0
      DO I=1,NBRICK
        INDCO = NODBRICK(:,I)
C  SET UP TETRAHEDRONS
        CALL MAKE_TETRA_48 (INDCO)

        ncell=NCELL+1
        if (nrtal.ne.nrad) ncltAL(ntact+1:ntet) = ncell

        IMATCH=0
        DO IQ=1,NQUAD

!  check side 1 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(11) == NODQUAD(9,IQ)) THEN
            INDF1 = (/ INDCO(21), INDCO(20), INDCO(19), 
     .                 INDCO(10), INDCO(1),  INDCO(2),
     .                 INDCO(3),  INDCO(12), INDCO(11) /)
            CALL BUBBLE (INDF1,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF1 == NSORQUAD)) THEN
              INMTIT(1,NTACT+1) = NQDSTS(IQ)
              INMTIT(1,NTACT+2) = NQDSTS(IQ)
              INMTIT(1,NTACT+3) = NQDSTS(IQ)
              INMTIT(1,NTACT+4) = NQDSTS(IQ)
              INMTIT(1,NTACT+5) = NQDSTS(IQ)
              INMTIT(1,NTACT+6) = NQDSTS(IQ)
              INMTIT(1,NTACT+7) = NQDSTS(IQ)
              INMTIT(1,NTACT+8) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+1,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+2,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+3,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+4,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+5,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+6,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+7,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+8,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 2 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(15) == NODQUAD(9,IQ)) THEN
            INDF2 = (/ INDCO(3),  INDCO(6),  INDCO(9), 
     .                 INDCO(18), INDCO(27), INDCO(24),
     .                 INDCO(21), INDCO(12), INDCO(15) /)
            CALL BUBBLE (INDF2,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF2 == NSORQUAD)) THEN
              INMTIT(1,NTACT+9) = NQDSTS(IQ)  
              INMTIT(1,NTACT+10) = NQDSTS(IQ) 
              INMTIT(1,NTACT+11) = NQDSTS(IQ)
              INMTIT(1,NTACT+12) = NQDSTS(IQ)
              INMTIT(1,NTACT+13) = NQDSTS(IQ)
              INMTIT(1,NTACT+14) = NQDSTS(IQ)
              INMTIT(1,NTACT+15) = NQDSTS(IQ)
              INMTIT(1,NTACT+16) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+9,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+10,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+11,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+12,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+13,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+14,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+15,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+16,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 3 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(17) == NODQUAD(9,IQ)) THEN
            INDF3 = (/ INDCO(9),  INDCO(8),  INDCO(7), 
     .                 INDCO(16), INDCO(25), INDCO(26),
     .                 INDCO(27), INDCO(18), INDCO(17) /)
            CALL BUBBLE (INDF3,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF3 == NSORQUAD)) THEN
              INMTIT(1,NTACT+17) = NQDSTS(IQ)
              INMTIT(1,NTACT+18) = NQDSTS(IQ)
              INMTIT(1,NTACT+19) = NQDSTS(IQ)
              INMTIT(1,NTACT+20) = NQDSTS(IQ)
              INMTIT(1,NTACT+21) = NQDSTS(IQ)
              INMTIT(1,NTACT+22) = NQDSTS(IQ)
              INMTIT(1,NTACT+23) = NQDSTS(IQ)
              INMTIT(1,NTACT+24) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+17,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+18,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+19,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+20,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+21,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+22,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+23,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+24,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 4 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(13) == NODQUAD(9,IQ)) THEN
            INDF4 = (/ INDCO(7),  INDCO(4),  INDCO(1), 
     .                 INDCO(10), INDCO(19), INDCO(22),
     .                 INDCO(25), INDCO(16), INDCO(13) /)
            CALL BUBBLE (INDF4,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF4 == NSORQUAD)) THEN
              INMTIT(1,NTACT+25) = NQDSTS(IQ) ! side 5
              INMTIT(1,NTACT+26) = NQDSTS(IQ) ! side 6
              INMTIT(1,NTACT+27) = NQDSTS(IQ)
              INMTIT(1,NTACT+28) = NQDSTS(IQ)
              INMTIT(1,NTACT+29) = NQDSTS(IQ)
              INMTIT(1,NTACT+30) = NQDSTS(IQ)
              INMTIT(1,NTACT+31) = NQDSTS(IQ)
              INMTIT(1,NTACT+32) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+25,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+26,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+27,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+28,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+29,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+30,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+31,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+32,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 5 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(5) == NODQUAD(9,IQ)) THEN
            INDF5 = (/ INDCO(3),  INDCO(2),  INDCO(1), 
     .                 INDCO(4),  INDCO(7),  INDCO(8),
     .                 INDCO(9),  INDCO(6),  INDCO(5) /)
            CALL BUBBLE (INDF5,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF5 == NSORQUAD)) THEN
              INMTIT(1,NTACT+33) = NQDSTS(IQ) 
              INMTIT(1,NTACT+34) = NQDSTS(IQ) 
              INMTIT(1,NTACT+35) = NQDSTS(IQ)
              INMTIT(1,NTACT+36) = NQDSTS(IQ)
              INMTIT(1,NTACT+37) = NQDSTS(IQ)
              INMTIT(1,NTACT+38) = NQDSTS(IQ)
              INMTIT(1,NTACT+39) = NQDSTS(IQ)
              INMTIT(1,NTACT+40) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+33,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+34,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+35,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+36,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+37,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+38,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+39,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+40,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 6 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(23) == NODQUAD(9,IQ)) THEN
            INDF6 = (/ INDCO(19), INDCO(20), INDCO(21), 
     .                 INDCO(24), INDCO(27), INDCO(26),
     .                 INDCO(25), INDCO(22), INDCO(23) /)
            CALL BUBBLE (INDF6,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF6 == NSORQUAD)) THEN
              INMTIT(1,NTACT+41) = NQDSTS(IQ)
              INMTIT(1,NTACT+42) = NQDSTS(IQ)
              INMTIT(1,NTACT+43) = NQDSTS(IQ)
              INMTIT(1,NTACT+44) = NQDSTS(IQ)
              INMTIT(1,NTACT+45) = NQDSTS(IQ)
              INMTIT(1,NTACT+46) = NQDSTS(IQ)
              INMTIT(1,NTACT+47) = NQDSTS(IQ)
              INMTIT(1,NTACT+48) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+41,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+42,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+43,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+44,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+45,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+46,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+47,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+48,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  more than 6 matches are impossible for a brick
!  so move to next brick
          IF (IMATCH == 6) EXIT
        END DO 

        ntact=ntet

      END DO

      WRITE (iunout,*) ' NUMBER OF SIDES OF TETRAHEDRON ON ',
     .            'NON-DEFAULT SURFACES ',COUNT(INMTIT>0)
      WRITE (iunout,*) ' NQUAD*8 ',NQUAD*8

      IF (.NOT.ALL(LUSED)) THEN
        WRITE (iunout,*) 
     .    ' SURFACES WITHOUT MATCHING BRICK FOUND IN IF0COP '
        DO IQ=1,NQUAD
          IF (.NOT.LUSED(IQ)) WRITE (iunout,*) IQ
        END DO
        CALL EXIT(1)
      END IF

!PB [1-sqrt(3/5)]-1,  0,  1- [....] 
      GK(1) = (1._DP - SQRT(3._DP/5._DP)) -1._DP
      GK(2) = 0._DP
      GK(3) = 1._DP - (1._DP - SQRT(3._DP/5._DP))

      DO I=1,1
        J = 0
        XVEC(1:27) = (/ (XTETRA(NODBRICK(K,I)), K=1,27) /)
        YVEC(1:27) = (/ (YTETRA(NODBRICK(K,I)), K=1,27) /)
        ZVEC(1:27) = (/ (ZTETRA(NODBRICK(K,I)), K=1,27) /)
        DO IR=1,3
           DO IS=1,3
              DO IT=1,3
                J=J+1
                CALL TRIQUAINT(TQVEC,GK(IR),GK(IS),GK(IT))
                XGAU=SUM(TQVEC*XVEC)*10._dp
                YGAU=SUM(TQVEC*YVEC)*10._dp
                ZGAU=SUM(TQVEC*ZVEC)*10._dp
                WRITE (iunout,'(2I6,3i3,3E20.10)') I,J,ir,is,it,
     .                XGAU,YGAU,ZGAU
              END DO
           END DO
        END DO
      END DO

!  CALCULATE PARTIAL DERIVATIVES DR/DX, DR/DY, DR/DZ, DS/DX, DS/DX, DS/DZ,
!  DT/DX, DT/DY AND DT/DZ FOR EACH BRICK (HEXAHEDRON) AT THE CENTER POINT

      ALLOCATE (DRSTDXYZ(3,3,NBRICK))
      CALL DFTRIQUA(DFVEC, 0._DP, 0._DP, 0._DP)

      DO I=1,NBRICK
        XVEC(1:27) = (/ (XTETRA(NODBRICK(K,I)), K=1,27) /)
        YVEC(1:27) = (/ (YTETRA(NODBRICK(K,I)), K=1,27) /)
        ZVEC(1:27) = (/ (ZTETRA(NODBRICK(K,I)), K=1,27) /)
        AJ(1,1) = SUM(XVEC*DFVEC(1,:))
        AJ(1,2) = SUM(YVEC*DFVEC(1,:))
        AJ(1,3) = SUM(ZVEC*DFVEC(1,:))
        AJ(2,1) = SUM(XVEC*DFVEC(2,:))
        AJ(2,2) = SUM(YVEC*DFVEC(2,:))
        AJ(2,3) = SUM(ZVEC*DFVEC(2,:))
        AJ(3,1) = SUM(XVEC*DFVEC(3,:))
        AJ(3,2) = SUM(YVEC*DFVEC(3,:))
        AJ(3,3) = SUM(ZVEC*DFVEC(3,:))

        AJM1T(1,1) = AJ(2,2)*AJ(3,3) - AJ(2,3)*AJ(3,2)
        AJM1T(1,2) = AJ(2,3)*AJ(3,1) - AJ(2,1)*AJ(3,3)
        AJM1T(1,3) = AJ(2,1)*AJ(3,2) - AJ(3,1)*AJ(2,2)
        AJM1T(2,1) = AJ(3,2)*AJ(1,3) - AJ(1,2)*AJ(3,3)
        AJM1T(2,2) = AJ(3,3)*AJ(1,1) - AJ(3,1)*AJ(1,3)
        AJM1T(2,3) = AJ(3,1)*AJ(1,2) - AJ(3,2)*AJ(1,1)
        AJM1T(3,1) = AJ(1,2)*AJ(2,3) - AJ(1,3)*AJ(2,2)
        AJM1T(3,2) = AJ(1,3)*AJ(2,1) - AJ(2,3)*AJ(1,1)
        AJM1T(3,3) = AJ(1,1)*AJ(2,2) - AJ(1,2)*AJ(2,1)

        DETI = 1._DP / (AJ(1,1)*AJ(2,2)*AJ(3,3) 
     .                + AJ(1,2)*AJ(2,3)*AJ(3,1)   
     .                + AJ(1,3)*AJ(2,1)*AJ(3,2)   
     .                - AJ(3,1)*AJ(2,2)*AJ(1,3)   
     .                - AJ(3,2)*AJ(2,3)*AJ(1,1)   
     .                - AJ(3,3)*AJ(2,1)*AJ(1,2))  

        AJM1 = TRANSPOSE(AJM1T) * DETI

        E = matmul(aj,ajm1)

        DRSTDXYZ(:,:,I) = AJM1
        if ((i==2333).or.(i==2340)) then
           write (iunout,*) ' i = ',i
           write (iunout,*) ' koordinaten '
           do j=1,27
              write (iunout,'(i6,4es12.4)') j,xvec(j),yvec(j),zvec(j)
           end do
           write (iunout,*) ' aj '
           write (iunout,'(3es12.4)') aj(1,1),aj(1,2),aj(1,3)
           write (iunout,'(3es12.4)') aj(2,1),aj(2,2),aj(2,3)
           write (iunout,'(3es12.4)') aj(3,1),aj(3,2),aj(3,3)
           write (iunout,*) ' ajm1 '
           write (iunout,'(3es12.4)') ajm1(1,1),ajm1(1,2),ajm1(1,3)
           write (iunout,'(3es12.4)') ajm1(2,1),ajm1(2,2),ajm1(2,3)
           write (iunout,'(3es12.4)') ajm1(3,1),ajm1(3,2),ajm1(3,3)
           write (iunout,*) ' determinate(aj) ',deti
!pb           if (i==2340) stop ' infcop '
        end if
      END DO


      DEALLOCATE (NODQUAD)
      DEALLOCATE (NQDSTS)
      DEALLOCATE (NQDSTP)
      DEALLOCATE (LUSED)

C
C
C  TRANSFER FLAGS
C
      NAINI=NAINB
C
C
      NLPLG=.FALSE.
      NLTET=.TRUE.
      LEVGEO=5
      NR1ST=NTET+1
      NLPOL=.FALSE.
      NP2ND=1
      NLTOR=.FALSE.
      NT3RD=1
      NR1TAL=NCLTAL(NTET)+1
      NP2TAL=NP2ND
      NT3TAL=NT3RD

      DO IN=NR1ST,NR1ST+NRADD
        NCLTAL(IN)=NCLTAL(NTET)+IN-NTET
      ENDDO

      CALL LEER(2)
      CALL HEADNG(' CASE REDEFINED IN COUPLE_TETRA: ',32)
      WRITE (iunout,*) 'NLPLG,NLFEM ',NLPLG,NLFEM
      WRITE (iunout,*) 'NLPOL       ',NLPOL
      WRITE (iunout,*) 'NR1ST,NP2ND ',NR1ST,NP2ND
      CALL LEER(2)
CTRIG E
C
      RETURN
C
C   GEOMETRY DEFINITION PART FINISHED
C
      ENTRY IF1COP
C
C   NOW READ THE PLASMA STATE GIVEN BY BRAAMS
C   AT PRESENT THE DATA COME FROM THE FILE FT31
C   FURTHERMORE: SCALING TO EIRENE UNITS AND INDEX MAPPING
C   STATEMENT NO. 2000 ---> 2999
C
C  IN CASE OF "SHORT CYCLE" THE PLASMA STATE IS TRANSFERRED VIA COMMON
C
      LSHORT=.FALSE.
      CALL LEER(1)
      WRITE (iunout,*) 'IF1COP CALLED '
      IF (NLPLAS) WRITE (iunout,*) 'PLASMA DATA EXPECTED ON BRAEIR'
      IF (.NOT.NLPLAS) 
     .  WRITE (iunout,*) 'PLASMA DATA EXPECTED ON FORT.31'
C  SKIP READING PLASMA, IF NLPLAS
      IF(NLPLAS) GOTO 2100
C
      GOTO 99991
C
C  IN CASE OF "SHORT CYCLE" OR TIME DEP. MODE
C  THE PLASMA STATE IS TRANSFERRED VIA COMMON
C  ONLY SCALING TO EIRENE UNITS AND INDEX MAPPING NEEDS TO BE DONE HERE
C
      ENTRY INTER1
      LSHORT=.TRUE.
      GOTO 2100
C
99991 CONTINUE
C
C
C  TRANSFER PROFILES
C
      IF (.NOT.(INDPRO(1).EQ.9.OR.INDPRO(2).EQ.9.OR.INDPRO(3).EQ.9.OR.
     .          INDPRO(4).EQ.9)) RETURN
C
C
      IF (NFLA.GT.NFL) THEN
        WRITE (iunout,*) ' PARAMETER ERROR DETECTED IN INFCOP '
        WRITE (iunout,*) ' NFLA MUST BE <= NFL'
        WRITE (iunout,*) ' NFLA,NFL = ',NFLA,NFL
        CALL EXIT(1)
      ENDIF

2100  CONTINUE
C
C  RESET 2D ARRAYS ONTO 1D EIRENE ARRAYS, RESCALE TO EIRENE UNITS
C  AND CONVERT BRAAMS VECTORS INTO CARTHESIAN EIRENE VECTORS
C
C  UNITS CONVERSION FACTORS
      T=1._DP/11600._DP
      V=1.
      VL=1.
CTRIG A
C  VACCUM DATA NEEDED FOR REGION OUTSIDE B2-MESH
      TVAC=0.02
      DVAC=1.E2_DP
      VVAC=0.
      BVAC=1.
CTRIG E
      DO 2105 IPLS=1,NPLSI
        D(IPLS)=FCTE(IPLS)*1.E-6_DP
        FL(IPLS)=FCTE(IPLS)
2105  CONTINUE
C
      ALLOCATE (NOSPEC(NPLSI))
      ALLOCATE (LUSED(NPLSI))

      ALLOCATE (DTEDX(NRTAL))
      ALLOCATE (DTEDY(NRTAL))
      ALLOCATE (DTEDZ(NRTAL))
      ALLOCATE (xlapla(NRTAL))

      OPEN (UNIT=31,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 31
      READ (31,'(A1000)') LINE

      ICOLUMN = 0
      IANF = 1
      NOSPEC = 0
      FORM='(A       )'
      LUSED=.FALSE.
      DO IPLS=1,NPLSI
        IF (LEN_TRIM(CDENMODEL(IPLS)) > 0) LUSED(IPLS)=.TRUE.
      END DO
      DO 
!  POSITION OF NEXT NON-BLANK CHARACTER
        I1 = VERIFY(LINE(IANF:),' ')
        IF (I1 == 0) EXIT
        I2 = INDEX(LINE(IANF+I1:),' ')
        ICOLUMN=ICOLUMN+1
        IF (ICOLUMN > 6) THEN
          IF (I2 < 10) THEN
            WRITE (FORM(3:3),'(I1)') I2
          ELSE
            WRITE (FORM(3:4),'(I2)') I2
          END IF
          SPECNAME=REPEAT(' ',72)
          READ(LINE(IANF+I1-1:IANF+I1+I2-2),FORM) SPECNAME(1:i2)
          SPECNAME=TRIM(SPECNAME)
          CALL UPPERCASE (SPECNAME)
          DO IPLS=1,NPLSI
            IF (SPECNAME == TEXTS(NSPAMI+IPLS)) THEN
              IF (NOSPEC(IPLS) == 0) THEN
                NOSPEC(IPLS) = ICOLUMN-6
                LUSED(IPLS)=.TRUE.
              ELSE
                WRITE (iunout,*) 
     .           ' ERROR! SPECIES ',SPECNAME,' FOUND TWICE '
              END IF
              EXIT
            END IF
          END DO
        END IF
        IANF=IANF+I1+I2-1
      END DO
      
      IF (.NOT. ALL(LUSED)) THEN
        WRITE (iunout,*) ' ERROR DETECTED IN INFCOP! '
        WRITE (iunout,*) 
     .    ' BULK PARTICLE DEFINITION IN EIRENE INPUT DOES NOT'
        WRITE (iunout,*) ' MATCH SPECIFICATIONS IN FIDAP PLASMA FILE '
        WRITE (iunout,*) ' CALCULATION ABANDONNED '
        CALL EXIT(1)
      END IF

      NREAD=maxval(NOSPEC(1:nplsi))
      ALLOCATE (TEF(NCOOR))
      ALLOCATE (DENF(0:MAX(NPLSI,NREAD),NCOOR))

      DO IC=1,NCOOR
        READ (31,'(A1000)') LINE
        READ (LINE,*) ICO,XCO,YCO,ZCO,TEF(IC),
     .                (DENF(IPLS,IC),IPLS=0,NREAD)
      END DO

      TEIN=TVAC
      TIIN=TVAC
      DIIN=DVAC
      DO ITET=1,NTET
        DIST1=1._DP/SQRT((XTETRA(NTECK(1,ITET))-XTCEN(ITET))**2+
     .                   (YTETRA(NTECK(1,ITET))-YTCEN(ITET))**2)
        DIST2=1._DP/SQRT((XTETRA(NTECK(2,ITET))-XTCEN(ITET))**2+
     .                   (YTETRA(NTECK(2,ITET))-YTCEN(ITET))**2)
        DIST3=1._DP/SQRT((XTETRA(NTECK(3,ITET))-XTCEN(ITET))**2+
     .                   (YTETRA(NTECK(3,ITET))-YTCEN(ITET))**2)
        DIST4=1._DP/SQRT((XTETRA(NTECK(4,ITET))-XTCEN(ITET))**2+
     .                   (YTETRA(NTECK(4,ITET))-YTCEN(ITET))**2)
        DIST=DIST1+DIST2+DIST3+DIST4
        TEIN(ITET) = T*(TEF(NTECK(1,ITET))*DIST1 + 
     .                  TEF(NTECK(2,ITET))*DIST2 +
     .                  TEF(NTECK(3,ITET))*DIST3 +
     .                  TEF(NTECK(4,ITET))*DIST4)/DIST
        DO IPLS=1,NPLSI
          IPL=NOSPEC(IPLS)
          IF (IPL == 0) CYCLE
          IPLSTI=MPLSTI(IPLS)
          TIIN(IPLSTI,ITET) = TEIN(ITET)
          DIIN(IPLS,ITET) = D(IPLS)*
     .                     (DENF(IPL,NTECK(1,ITET))*DIST1 + 
     .                      DENF(IPL,NTECK(2,ITET))*DIST2 +
     .                      DENF(IPL,NTECK(3,ITET))*DIST3 +
     .                      DENF(IPL,NTECK(4,ITET))*DIST4)/DIST
        END DO
      END DO

      write (iunout,*) ' temax = ',maxval(tein(1:NTET))
      write (iunout,*) ' temin = ',minval(tein(1:NTET))
!  koordinaten der gaussknoten in r,s,t
!PB [1-sqrt(3/5)]-1,  0,  1- [....] 
      GK(1) = (1._DP - SQRT(3._DP/5._DP)) -1._DP
      GK(2) = 0._DP
      GK(3) = 1._DP - (1._DP - SQRT(3._DP/5._DP))

        

      dtedx = 0._dp
      dtedy = 0._dp
      dtedz = 0._dp
      xlapla = 0._dp
!   eckpunkte des bricks (27 pro brick)
      do ir=1,nr1tal-1
        XVEC(1:27) = (/ (XTETRA(NODBRICK(K,IR)), K=1,27) /)
        YVEC(1:27) = (/ (YTETRA(NODBRICK(K,IR)), K=1,27) /)
        ZVEC(1:27) = (/ (ZTETRA(NODBRICK(K,IR)), K=1,27) /)
        TEVEC(1:27) = T*(/ (TEF(NODBRICK(K,IR)), K=1,27) /)
!   calculate gradient T in each cell (centre of brick)
        call gradf (0._dp, 0._dp, 0._dp, xvec, yvec, zvec, tevec,
     .              dtedx(ir), dtedy(ir), dtedz(ir))
!   gradient length in each cell, for comparison with photon mean free path
        TEDTEDX(IR) =  TEVEC(14) / DTEDX(IR)
        TEDTEDY(IR) =  TEVEC(14) / DTEDY(IR)
        TEDTEDZ(IR) =  TEVEC(14) / DTEDZ(IR)

!   calculate laplace T in each cell (centre of brick)
        xlapla(ir) = flaplace(0._dp, 0._dp, 0._dp, 
     .                         xvec, yvec, zvec, tevec)
      end do

      DEALLOCATE (TEF)
      DEALLOCATE (DENF)
      DEALLOCATE (NOSPEC)
      DEALLOCATE (LUSED)

C
C
      RETURN
C
2999  CONTINUE
C
C  PLASMA PROFILES ARE NOW READ IN
C
      ENTRY IF2COP(ITARG)
      IF (ITARG.GT.NTARGI) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'SOURCE DATA FOR STRATUM ISTRA= ',ITARG
        WRITE (iunout,*) 
     .    'CANNOT BE DEFINED IN IF2COP. CHANGE INDSRC(ISTRA)'
        CALL LEER(1)
        RETURN
      ENDIF
C
C  NEXT DEFINE FLUXES, TEMPERATURES AND VELOCITIES AT THE TARGETS
C  (FLUXES IN AMP/(CM ALONG TARGET), TEMPERATURES IN EV, VELOCITIES IN CM/SEC)
C   FNIXB*FL (FNIYB*FL) ARE GIVEN IN AMP
C  STATEMENT NO 3000 ---> 3999
C
3000  CONTINUE
C
      
      TESTEP(ITARG,:)=0.D0
      TISTEP(:,ITARG,:)=0.D0
      DISTEP(:,ITARG,:)=0.D0
      VXSTEP(:,ITARG,:)=0.D0
      VYSTEP(:,ITARG,:)=0.D0
      VZSTEP(:,ITARG,:)=0.D0
      FLSTEP(:,ITARG,:)=0.D0
      DO IG=1,NRWL(ITARG)
        IN=IRSTEP(ITARG,IG)
        TESTEP(ITARG,IG)=TEIN(IN)
        TISTEP(1:NPLSTI,ITARG,IG)=TIIN(1:NPLSTI,IN)
        DISTEP(1:NPLSI,ITARG,IG)=DIIN(1:NPLSI,IN)
        VXSTEP(1:NPLSV,ITARG,IG)=VXIN(1:NPLSV,IN)
        VYSTEP(1:NPLSV,ITARG,IG)=VYIN(1:NPLSV,IN)
        VZSTEP(1:NPLSV,ITARG,IG)=VZIN(1:NPLSV,IN)
        DO IPLS=1,NPLSI
          IPLSTI = MPLSTI(IPLS)
          CSP(IPLS)=CVEL2A*SQRT((TIIN(IPLSTI,IN)+TEIN(IN))/
     .              RMASSP(IPLS))
          FLSTEP(IPLS,ITARG,IG)=ELCHA*0.5*DIIN(IPLS,IN)*CSP(IPLS)
        END DO
      END DO
      nrwl(itarg)=nrwl(itarg)+1
C
      IF (TRCSOU) CALL LEER(2)
C
C  INITALIZE FUNCTION STEP (FOR RANDOM SAMPLING ALONG TARGET)
C  SET SOME SOURCE PARAMETERS EXPLICITLY TO ENFORCE INPUT CONSISTENCY
C
      IIPLS=1
      IEPLS=NPLSI
      FLUX(ITARG)=STEP(IIPLS,IEPLS,NRWL(ITARG),ITARG)
C
      NLPLS(ITARG)=.TRUE.
      NLATM(ITARG)=.FALSE.
      NLMOL(ITARG)=.FALSE.
      NLION(ITARG)=.FALSE.
C
      NLSRF(ITARG)=.TRUE.
      NLPNT(ITARG)=.FALSE.
      NLLNE(ITARG)=.FALSE.
      NLVOL(ITARG)=.FALSE.
      NLCNS(ITARG)=.FALSE.
C
      NSRFSI(ITARG)=1
      INDIM(1,ITARG)=4
      I4=IDEZ(INT(SORLIM(1,ITARG)),4,4)
      SORLIM(1,ITARG)=I4*1000+104
      SORIND(1,ITARG)=ITARG
C  IN CASE INDIM=4: INSOR,INDGRD... ARE REDUNDANT
      NRSOR(1,ITARG)=-1
      NPSOR(1,ITARG)=-1
      IF (INDSRC(ITARG).LT.6) THEN
        WRITE (iunout,*) 'MESSAGE FROM IF2COP: '
        WRITE (iunout,*) 'SOURCE STRENGTH AND SPATIAL DISTRIBUTION FOR '
        WRITE (iunout,*) 'STRATUM ',ISTRA,' MODIFIED.'
        CALL MASR1('FLUX=   ',FLUX(ISTRA))
        WRITE (iunout,*) 'USE STEP FUNCTION ISTEP= ',ITARG,
     .    ' FROM BLOCK 14'
      ENDIF
C
      IF (INDSRC(ITARG).EQ.6) THEN
C  DEFINE SOURCE FOR TARGET RECYCLING STRATUM ITARG
C  ASSUME NOW: ITARG=ISTRA
C  DEFAULTS ARE ALREADY SET IN SUBR. INPUT.
C
        CALL FTCRI(ITARG,CITARG)
        TXTSOU(ITARG)= 'SURFACE RECYCLING SOURCE NO.'//CITARG
        NPTS(ITARG)=NPTC(ITARG,1)
        NINITL(ITARG)=ITARG*1001
        NSPEZ(ITARG)=-1
        SORIFL(1,ITARG)=NIFLG(ITARG,1)
        SORWGT(1,ITARG)=1.
        IF (NIXY(ITARG,1).EQ.1) THEN
C TARGET RECYCLING SOURCE AT POLOIDAL SURFACE NPES
          NEMODS(ITARG)=3
          NAMODS(ITARG)=1
          SORENI(ITARG)=3.
          SORENE(ITARG)=0.5
        ELSEIF (NIXY(ITARG,1).EQ.2) THEN
C WALL RECYCLING SOURCE AT RADIAL SURFACE NPES
          NEMODS(ITARG)=2
          NAMODS(ITARG)=1
          SORENI(ITARG)=2.
          SORENE(ITARG)=0.
        ENDIF
C
C  SORAD1,...: USE POLYGON MESH, IE. SORAD1,... ARE REDUNDANT.
C
C  VELOCITY SPACE DISTRIBUTION
        SORCOS(ITARG)=1.
        SORMAX(ITARG)=0.
C
C
C  DO 2028 LOOP FROM SUBR. INPUT
        THMAX=MAX(0._DP,MIN(PIHA,SORMAX(ITARG)*DEGRAD))
        IF (NAMODS(ITARG).EQ.1) THEN
          RP1=SORCOS(ITARG)+1.
          SORCOS(ITARG)=1./RP1
          IF (ABS(COS(THMAX)).LE.EPS10) THEN
            SORMAX(ITARG)=1.
          ELSE
            SORMAX(ITARG)=1.-COS(THMAX)**RP1
          ENDIF
        ELSEIF (NAMODS(ITARG).EQ.2) THEN
          SORCOS(ITARG)=SORCOS(ITARG)*DEGRAD
          SORMAX(ITARG)=THMAX
        ENDIF
        NLSYMT(0)=NLSYMT(0).AND.NLSYMT(ITARG)
        NLSYMP(0)=NLSYMP(0).AND.NLSYMP(ITARG)
C
      ENDIF
C
C  SOURCE DEFINITION FOR TARGET RECYCLING STRATUM ITARG COMPLETED
C
3999  CONTINUE
C
C  TARGET DATA ARE DEFINED NOW
C
      DO 5000 IG=1,NGITT
        ELSTEP(:,ITARG,IG)=0.
5000  CONTINUE
C
C  COMPUTE EXACT SURFACE ENERGY FLUXES FOR COMPARISON WITH SAMPLED
C  E-FLUX "ETOTP". THIS IS ONLY FOR DIAGNOSTICS PURPOSES
C  E.G. TO CHECK CONSISTENCY OF BOUNDARY CONDITIONS
C  STATEMENT NO. 6000 ---> 6500
C
      IF (.NOT.TRCSOU) GOTO 6500
C
      write (iunout,*) ' exact surface energy flux computation '
      write (iunout,*) ' not yet available of tetrahedrons '
      GOTO 6500

 4812 continue
      EEMAX=0.
      EESHT=0.
C
      DO 6011 IG=1,NRWL(ITARG)-1
        OR=ORI(ITARG,IG)
C  COMPUTE SHEATH POTENTIAL ESHT(ITARG,IG)
C  USE ALL NPLSI SPECIES, NOT JUST IFL=NSPZI,NSPZE
        ESHT(ITARG,IG)=0.D0
        IF (NEMODS(ITARG).EQ.3.OR.NEMODS(ITARG).EQ.5.OR.
     .      NEMODS(ITARG).EQ.7) THEN
          IF (IGSTEP(ITARG,IG).GT.200000) THEN
            ITRI=IRSTEP(ITARG,IG)
            NPES=IGSTEP(ITARG,IG)-200000
            DO 6005 IPL=1,NPLSI
              IPLV=MPLSV(IPL)
              PM1=(PTRIX(NPES,ITRI)*VXSTEP(IPLV,ITARG,IG)+
     .             PTRIY(NPES,ITRI)*VYSTEP(IPLV,ITARG,IG))*OR
              VPZ=VZSTEP(IPLV,ITARG,IG)
              VP(IPL)=SQRT(PM1**2+VPZ**2)
              DI(IPL)=DISTEP(IPL,ITARG,IG)
6005        CONTINUE
            TE=TESTEP(ITARG,IG)
            CUR=0.
            GAMMA=0.
            ESHT(ITARG,IG)=SHEATH(TE,DI,VP,NCHRGP,GAMMA,CUR,NPLSI,
     .                           -ITARG)
          ELSEIF (IGSTEP(ITARG,IG).LT.200000) THEN
            ITRI=IRSTEP(ITARG,IG)
            NPES=IGSTEP(ITARG,IG)-100000
            ESHT(ITARG,IG)= 0
          ENDIF
        ENDIF
C
        IF (IGSTEP(ITARG,IG).LT.200000) GOTO 6010
        ITRI=IRSTEP(ITARG,IG)
        NPES=IGSTEP(ITARG,IG)-200000
        DO 6009 IPLS=1,NPLSI
          IF (FLSTEP(IPLS,ITARG,IG).EQ.0.D0) GOTO 6009
          IPLSTI=MPLSTI(IPLS)
          IPLSV=MPLSV(IPLS)
          VT=SQRT(2.*TISTEP(IPLSTI,ITARG,IG)/BMASS(IPLS))*CVEL2A
C  VELOCITY COMPONENT NORMAL TO TARGET SURFACE
C  I.E., POLOIDAL COMPONENT V-POL
C  ASSUMING ORTHOGONAL TARGET
          PM1=(PTRIX(NPES,ITRI)*VXSTEP(IPLSV,ITARG,IG)+
     .         PTRIY(NPES,ITRI)*VYSTEP(IPLSV,ITARG,IG))*OR
C  VELOCITY COMPONENT PARALLEL TO TARGET SURFACE
C  I.E., RADIAL PLUS TOROIDAL COMPONENT, V-RAD + V-TOR
C  AGAIN: ASSUMING ORTHOGONAL TARGET
          VPX=VXSTEP(IPLSV,ITARG,IG)-PM1*PPLNX(IY,NPES)*OR
          VPY=VYSTEP(IPLSV,ITARG,IG)-PM1*PPLNY(IY,NPES)*OR
          VPZ=VZSTEP(IPLSV,ITARG,IG)-0.
          PN1=SQRT(VPX**2+VPY**2+VPZ**2)
          PERW=0.
          PARW=0.
          IF (VT.GT.0.) THEN
            PERW=PM1/VT
            PARW=PN1/VT
          ENDIF
C
          CS=SQRT((1.*TISTEP(IPLSTI,ITARG,IG)+
     .                TESTEP(ITARG,IG))/BMASS(IPLS))*CVEL2A
C THE MACH NUMBER BOUNDARY CONDITION ONLY AFFECTS THE PARALLEL TO B
C MOMENTUM, I.E., NOT THE RADIAL VELOCITY
          VTEST=SQRT(PM1**2+VPZ**2)
          VTEST=VTEST/(CS+EPS60)
          VR=SQRT(VPX**2+VPY**2)
          WRITE (iunout,*) 'IPLS,ITARG,IG,MACH ',IPLS,ITARG,IG,VTEST
C         WRITE (iunout,*) 'POL., TOR., RAD. ',PM1,VPZ,VR
          CALL LEER(1)
C TARGET ENERGY FLUXES
          DRR=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
          IF (NEMODS(ITARG).EQ.1) THEN
            EADD=SORENI(ITARG)
          ELSEIF (NEMODS(ITARG).EQ.2.OR.NEMODS(ITARG).EQ.3) THEN
            EADD=SORENI(ITARG)*TISTEP(IPLSTI,ITARG,IG)+SORENE(ITARG)*
     .             TESTEP(ITARG,IG)
          ELSEIF (NEMODS(ITARG).GE.4) THEN
            PERWI=PERW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
            PARWI=PARW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
            EADD=EMAXW(TISTEP(IPLSTI,ITARG,IG),PERWI,PARWI)
          ENDIF
          ESUM=EADD*FLSTEP(IPLS,ITARG,IG)
          ELSTEP(IPLS,ITARG,IG)=ELSTEP(IPLS,ITARG,IG)+ESUM
          EEMAX=EEMAX+ESUM*DRR
C  ADD SHEATH ACCELERATION
          EADD=NCHRGP(IPLS)*ESHT(ITARG,IG)
          ESUM=EADD*FLSTEP(IPLS,ITARG,IG)
          EESHT=EESHT+ESUM*DRR
          ELSTEP(IPLS,ITARG,IG)=ELSTEP(IPLS,ITARG,IG)+ESUM
6009    CONTINUE
        GOTO 6011
6010    CONTINUE
C  TO BE WRITTEN
6011  CONTINUE
C
      CALL LEER(1)
      WRITE (iunout,*) 'TARGET DATA: TARGET NO. ITARG=ISTRA= ',ITARG
      WRITE (iunout,*) 'IG, ARC, P-FLUX, E-FLUX, TE, TI, SHEATH/TE'
      DO 6100 IG=1,NRWL(ITARG)-1
        WRITE (iunout,'(1X,I4,1P,6E11.3)')
     .             IG,RRSTEP(ITARG,IG),FLSTEP(0,ITARG,IG),
     .             ELSTEP(0,ITARG,IG),
     .             TESTEP(ITARG,IG),TISTEP(1,ITARG,IG),
     .             ESHT(ITARG,IG)/(TESTEP(ITARG,IG)+EPS60)
6100  CONTINUE
      WRITE (iunout,'(1X,I4,1P,1E11.3)') NRWL(ITARG),
     .                                 RRSTEP(ITARG,NRWL(ITARG))
      CALL MASR1 ('EEMAX    ',EEMAX)
      CALL MASR1 ('EESHT    ',EESHT)
C
      ETOT=EEMAX+EESHT
      EFLX(ITARG)=EEMAX+EESHT
      WRITE (iunout,*) 'PARTICLE FLUX(IPLS), IPLS=1,NPLSI '
      WRITE (iunout,'(1X,1P,6E12.4)') (FLTOT(ISP,ITARG),ISP=1,NPLSI)
      CALL LEER(1)
      WRITE (iunout,*) 'ENERGY FLUX '
      WRITE (iunout,'(1X,1P,1E12.4)') EFLX(ITARG)
      CALL LEER(2)
C
6300  CONTINUE
C
C
C  SET SOME OTHER DATA SPECIFIC FOR EIRENE CODE REQUIREMENTS
C  STATEMENT NO. 6500 ---> 6999
C
6500  CONTINUE
C
C
      RETURN
999   CONTINUE
      WRITE (iunout,*) 'ERROR IN IF2COP: NGITT TOO SMALL '
      CALL EXIT(1)
      RETURN
C
C
      ENTRY IF3COP(ISTRAA,ISTRAE,NEW_ITER)
C
C
      WRITE (iunout,*) ' IF3COP IS CALLED, ISTRAA,ISTRAE '
      WRITE (iunout,*) ISTRAA,ISTRAE
      LSHORT=.FALSE.
      LSTOP=.TRUE.
      IFIRST=0
      GOTO 99992
C
      ENTRY INTER3(LSTP,IFRST,ISTRAA,ISTRAE,NEW_ITER)
C
C  ENTRY FOR SHORT CYCLE FROM SUBR. EIRSRT
C
C  IFIRST=0: RESTORE DATA FROM A PREVIOUS EIRENE RUN, SET REFERENCE
C            DATA FOR "STOP-CRITERION" SNIS,SEES,SEIS
C  IFIRST>0: MODIFY SOURCE TERMS ACCORDING TO NEW PLASMA CONDITIONS,
C            COMPARE INTEGRALS WITH SNIS,...., AND DECIDE TO STOP OR
C            CONTINUE SHORT CYCLE (LSTOP)
C
      LSHORT=.TRUE.
      LSTOP=LSTP
      IFIRST=IFRST
C
99992 CONTINUE
C
      DO 10000 ISTRAI=ISTRAA,ISTRAE
C
        IF (XMCP(ISTRAI).LE.1.) GOTO 10000
C
        IF (LSHORT) GOTO 7000
C
C  READ DATA FROM STRATUM NO. ISTRAI BACK INTO WORKING SPACE
C  IF REQUIRED
C
        IF (ISTRAI.EQ.IESTR) THEN
C  NOTHING TO BE DONE
        ELSEIF ((NFILEN.EQ.1.OR.NFILEN.EQ.2).AND.IESTR.NE.ISTRAI) THEN
          IESTR=ISTRAI
          CALL RSTRT(ISTRAI,NSTRAI,NESTM1,NESTM2,NADSPC,
     .               ESTIMV,ESTIMS,ESTIML,
     .               NSDVI1,SDVI1,NSDVI2,SDVI2,
     .               NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .               NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .               NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .               NSIGI_SPC,TRCFLE)
        ELSE
          WRITE (iunout,*) 'ERROR IN INFCOP: STRATUM ISTRAI= ',ISTRAI
          WRITE (iunout,*) 'IS NOT AVAILABLE. EXIT CALLED'
          CALL EXIT(1)
        ENDIF
C
C  DATA TRANSFER BACK FROM EIRENE TO EXTERNAL CODE
C  STATEMENT NO 7000 ---> 7999
C
7000    CONTINUE
C
C  SCALE SURFACE SOURCES PER UNIT FLUX, FOR OTHER SOURCES USE
C  EIRENE SCALINGS
        IF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).NE.0.) THEN
C  FLUX FROM EIRENE TO PLASMA CODE: NEGATIVE
          FLX=-WTOTP(0,ISTRAI)
          FLXI=1./FLX
        ELSEIF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).EQ.0.) THEN
          WRITE (iunout,*) 'NO FLUX FROM STRATUM NO. ISTRAI= ',ISTRAI
          WRITE (iunout,*) 'NO DATA RETURNED FOR THIS STRATUM'
          GOTO 7999
        ELSEIF (ISTRAI.GT.NTARGI) THEN
          FLXI=1.
        ENDIF
C
        IF (.NOT.LSHORT) GOTO 7400
C  SHORT LOOP CORRECTION FINISHED
C
7400    CONTINUE
C
C
C  ADD CONTRIBUTIONS FROM VOLUME RECOMBINATION SOURCE
C
C
        EPEL = 0.D0
        EMIS = 0._DP
        ABSORB = 0._DP
        VOLSUM = 0._DP
        ABSAMP = 0._DP
        AB1SMP = 0._DP
        AB2SMP = 0._DP
        EMSAMP = 0._DP
        EM1SMP = 0._DP
        EM2SMP = 0._DP
        RIX = 0._DP
        RIY = 0._DP
        RIZ = 0._DP
        XKAP = 0._DP
        YKAP = 0._DP
        ZKAP = 0._DP

        IF (NLVOL(ISTRAI)) THEN
C
C  EPEL: CALCULATED EMISSION (NOT SAMPLED)
C
          RECTOT = 0._DP
          DO ISRFSI=1,NSRFSI(ISTRAI)
            ISR=ISRFSI
            ISTEP=SORIND(ISR,ISTRAI)
            DO IPLS=1,NPLSI
              DO IIRC=1,NPRCI(IPLS)
              IRRC=LGPRC(IPLS,IIRC)
                IF ((ISTEP.EQ.0).OR.(ISTEP.EQ.IRRC)) THEN
                  SUMN=0.0
                  SUMEE=0.0
                  DO ITET=1,NTET
                    INC=NCLTAL(ITET)
                    IF (NSTORDR >= NRAD) THEN
                      RECADD=-TABRC1(IRRC,ITET)*DIIN(IPLS,ITET)*ELCHA
                      EEADD=  EELRC1(IRRC,ITET)*DIIN(IPLS,ITET)*ELCHA
                    ELSE
                      RECADD=-FTABRC1(IRRC,ITET)*DIIN(IPLS,ITET)*ELCHA
                      EEADD=  FEELRC1(IRRC,ITET)*DIIN(IPLS,ITET)*ELCHA
                    END IF
                    SUMN=SUMN+RECADD*VOL(ITET)
                    EPEL(INC)=EPEL(INC)+EEADD*VOL(ITET)
                    SUMEE=SUMEE+EEADD*VOL(ITET)
                  END DO         ! ITET
                  RECTOT = RECTOT + SUMN
                  WRITE (iunout,*) 'IPLS,IRRC ',IPLS,IRRC
                  CALL MASR2('SUMN, SUMEE     ',
     .                        SUMN,SUMEE)
                END IF
              END DO             ! IIRC
            END DO               ! IPLS 
          END DO                 ! ISRFSI
          do i=1,nsbox_tal
             epel(i) = epel(i)/voltal(i)
          END DO                 
        END IF

C  NEW VERSION FOR KAPPA_RAD
C  KAPPA FROM  NET SOURCE, THICK FRACTION, AND LAPLACE T 
C  (ANY T WILL DO, ALL ARE THE SAME)

        do i=1,nsbox_tal
           COPV(10,I)=-(COPV(3,I)+copv(4,I))/(xlapla(I)+EPS60)
           COPV(11,I)=xlapla(i)

! if kappa is negative put sources for this cell into thin fraction 
! note:  the variance for the total (sigma_cop(nspvi+1)) is not affected by this
!        but the individual variances sigma_cop(3...6) are not correct anymore.
            LSHIFT=COPV(10,I) < 0._DP
            IF (LSHIFT) THEN
              COPV(5,I)  = COPV(5,I) + COPV(3,I)
              COPV(3,I)  = 0._DP
              COPV(6,I)  = COPV(6,I) + COPV(4,I)
              COPV(4,I)  = 0._DP
!   set radiative vector flux I  of thick fraction to zero in this cell
              COPV(7,I)  = 0._DP
              COPV(8,I)  = 0._DP
              COPV(9,I)  = 0._DP
!   set kappa to zero in this cell
              COPV(10,I) = 0._DP
            END IF
        end do


        do i=3,13
          HELPAR(1:NSBOX_TAL) = COPV(I,1:NSBOX_TAL)
          CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
     .                 COPVI(I,ISTRAI),
     .                 NR1TAL,NP2TAL,NT3TAL,NBMLT)
          COPV(I,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)
        end do


! NEXT:  FILL ADDV, FOR PLOTTING, VARIANCES AND DIAGNOSTICS ARRAYS

!  addv(1):  calculated emission
!  addv(2):  calculated emission - sampled absorption
!  addv(3):  standard deviation (%) for addv(2)
!  addv(4):  sampled emission - sampled absorption
!  addv(5):  standard deviation (%) for sampled emission - sampled absorption
!! NEXT: INTERMEDIATE DIAGNOSTICS,  6-8: CURRENTLY OUT
!!  addv(6):  dT/dx,  and later: div(vec(I)
!!  addv(7):  dT/dy
!!  addv(8):  dT/dz

        ADDV(1,:) = EPEL
!       ADDV(2,:) = EPEL + EPHPHT = EPEL + COPV(2)
        ADDV(2,:) = EPEL + COPV(2,:)

!  total emission is calculated. Standard deviation in absolute units
!  rescale sigma(,..) 1. )to absolute values, then
!                     2.) from total absorption to net emission.
!                     since total emission is an additive constant,
!                     step 2.) is not necessary
        WHERE (ADDV(2,:) /= 0._DP)
          ADDV(3,:) = SIGMA_COP(2,:)/100._DP*ABS(COPV(2,:)) 
        ELSE WHERE
          ADDV(3,:) = 0._DP
        END WHERE
!

!  total emission is sampled. 
        ADDV(4,:) = COPV(3,:) + COPV(4,:) + COPV(5,:) + COPV(6,:)
!  rescale sigma(,..) to absolute values 
        ADDV(5,:) = SIGMA_COP(NCPVI+1,:)/100._DP*ABS(ADDV(4,:))
!
        HELPAR(1:NSBOX_TAL) = ADDV(1,1:NSBOX_TAL)
        CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
     .               ADDVI(1,ISTRAI),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        ADDV(1,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)
        
        HELPAR(1:NSBOX_TAL) = ADDV(2,1:NSBOX_TAL)
        CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
     .               ADDVI(2,ISTRAI),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        ADDV(2,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)

!  STD DEVIATION OF CELL AVERAGE, ABSOLUTE VALUES
        ADDVI(3,ISTRAI) = SGMS_COP(2)/100._DP*ABS(COPVI(2,ISTRAI))
        
        HELPAR(1:NSBOX_TAL) = ADDV(4,1:NSBOX_TAL)
        CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
     .               ADDVI(4,ISTRAI),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        ADDV(4,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)

!  STD DEVIATION OF CELL AVERAGE, ABSOLUTE VALUES
        ADDVI(5,ISTRAI) = SGMS_COP(NCPVI+1)/100._DP*
     .                    ABS(ADDVI(4,ISTRAI))

!       ADDV(6,1:NSBOX_TAL) = DTEDX(1:NSBOX_TAL)
!       ADDV(7,1:NSBOX_TAL) = DTEDY(1:NSBOX_TAL)
!       ADDV(8,1:NSBOX_TAL) = DTEDZ(1:NSBOX_TAL)

!       HELPAR(1:NSBOX_TAL) = ADDV(6,1:NSBOX_TAL)
!       CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
!    .               ADDVI(6,ISTRAI),
!    .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
!       ADDV(6,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)
!       HELPAR(1:NSBOX_TAL) = ADDV(7,1:NSBOX_TAL)
!       CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
!    .               ADDVI(7,ISTRAI),
!    .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
!       ADDV(7,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)
!       HELPAR(1:NSBOX_TAL) = ADDV(8,1:NSBOX_TAL)
!       CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
!    .               ADDVI(8,ISTRAI),
!    .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
!       ADDV(8,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)

C
        IF (.NOT.LSYMET) GOTO 7500
C
C  SECONDLY SYMMETRISE EIRENE ARRAYS ACCORDING TO SYMMETRY IN MODEL
C
C
C   THIRDLY WRITE EIRENE ARRAYS (1D) ONTO BRAAMS ARRAYS (2D)
C   AND RESCALE TO PROPER UNITS: #/CELL/STRATUM FLUX
C   # STANDS FOR PARTICLES (SNI), MOMENTUM (SMO)
C   AND ENERGY (SEE,SEI)
C
7500    CONTINUE
C
C   NEXT:
C   IF LSHORT: CRITERION TO STOP SHORT CYCLE,
C   IF NOT LSHORT: RESCALE SURFACE SOURCE STRATA
C                  UNITS: # PER UNIT TARGET PLATE FLUX
C
C
C
C   THIRDLY:
C   INDEX MAPPING BACK TO BRAAMS IMPLEMENTATION OF LINDA GEOMETRY
C
C
7700    CONTINUE
C
7999    CONTINUE

c  extrapolate tallies from cell averages to fidap-vertices (27 for each brick)

        DO ITET=1,NTET
          INC=NCLTAL(ITET)
          DO J=1,4
            ipo=nteck(j,itet)
            dist1=1._dp/sqrt((xtetra(ipo)-xtcen(itet))**2 +
     .                       (ytetra(ipo)-ytcen(itet))**2 +
     .                       (ztetra(ipo)-ztcen(itet))**2) 
            volsum(ipo)=volsum(ipo)+dist1
            EMIS(ipo) = EMIS(ipo) + EPEL(NCLTAL(ITET))*dist1
            ABSORB(ipo) = ABSORB(ipo)+
     .                    EPHPHT(NCLTAL(ITET))*dist1
            EMSAMP(ipo) = EMSAMP(IPO) + COPV(1,NCLTAL(ITET))*DIST1
            ABSAMP(ipo) = ABSAMP(IPO) + COPV(2,NCLTAL(ITET))*DIST1
            EM1SMP(ipo) = EM1SMP(IPO) + COPV(3,NCLTAL(ITET))*DIST1
            AB1SMP(ipo) = AB1SMP(IPO) + COPV(4,NCLTAL(ITET))*DIST1
            EM2SMP(ipo) = EM2SMP(IPO) + COPV(5,NCLTAL(ITET))*DIST1
            AB2SMP(ipo) = AB2SMP(IPO) + COPV(6,NCLTAL(ITET))*DIST1
            RIX(ipo)    = RIX(IPO) + COPV(7,NCLTAL(ITET))*DIST1
            RIY(ipo)    = RIY(IPO) + COPV(8,NCLTAL(ITET))*DIST1
            RIZ(ipo)    = RIZ(IPO) + COPV(9,NCLTAL(ITET))*DIST1
!   kappa_rad from grad Te:  out
!   copv(11,11,12) have now another meaning
!pb         xkap(ipo)   = xkap(IPO) + COPV(10,NCLTAL(ITET))*DIST1
!pb         ykap(ipo)   = ykap(IPO) + COPV(11,NCLTAL(ITET))*DIST1
!pb         zkap(ipo)   = zkap(IPO) + COPV(12,NCLTAL(ITET))*DIST1
!   kappa_rad from laplace Te
            xkap(ipo)   = xkap(IPO) + COPV(10,NCLTAL(ITET))*DIST1
          END DO                 ! J
        END DO                   ! ITET

        where (abs(volsum(1:ncoor)) > 1.D-10) 
          EMIS(1:ncoor) = EMIS(1:ncoor)/volsum(1:ncoor)
          ABSORB(1:ncoor) = ABSORB(1:ncoor)/volsum(1:ncoor)
          EMSAMP(1:ncoor) = EMSAMP(1:ncoor)/volsum(1:ncoor)
          ABSAMP(1:ncoor) = ABSAMP(1:ncoor)/volsum(1:ncoor)
          EM1SMP(1:ncoor) = EM1SMP(1:ncoor)/volsum(1:ncoor)
          EM2SMP(1:ncoor) = EM2SMP(1:ncoor)/volsum(1:ncoor)
          AB1SMP(1:ncoor) = AB1SMP(1:ncoor)/volsum(1:ncoor)
          AB2SMP(1:ncoor) = AB2SMP(1:ncoor)/volsum(1:ncoor)
          RIX(1:ncoor) = RIX(1:ncoor)/volsum(1:ncoor)*elcha
          RIY(1:ncoor) = RIY(1:ncoor)/volsum(1:ncoor)*elcha
          RIZ(1:ncoor) = RIZ(1:ncoor)/volsum(1:ncoor)*elcha
!pb       XKAP(1:ncoor) = XKAP(1:ncoor)/volsum(1:ncoor)
!pb       YKAP(1:ncoor) = YKAP(1:ncoor)/volsum(1:ncoor)
!pb       ZKAP(1:ncoor) = ZKAP(1:ncoor)/volsum(1:ncoor)
          XKAP(1:ncoor) = XKAP(1:ncoor)/volsum(1:ncoor)
        end where

! TALLIES FOR FIDAP INTERFACE: DONE
!  write output file for FIDAP, stream 35


        IF (ISTRAI == 1) THEN
          WRITE (35,*) NSTRAI
          WRITE (35,*) NCOORD
        END IF
        WRITE (35,*) TXTSOU(ISTRAI)
        WRITE (35,'(A9,1X,10A20)') 'NO. POINT','EMIS_1','ABSORPTION_1',
     .                            'EMIS_2','ABSORBTION_2',
     .                            'Kappa','I_x','I_y','I_z'
        DO IR=1,NCOORD
          WRITE (35,'(I9,1X,10ES20.10)') IR, EM1SMP(IR),AB1SMP(IR),
     .                                  EM2SMP(IR),AB2SMP(IR),
     .                                  xkap(ir),
     .                                  RIX(IR),RIY(IR),RIZ(IR)
        END DO
        WRITE (35,*)
C
C  DATA TRANSFER BACK TO PLASMA CODE FINISHED FOR STRATUM NO. ISTRAI
C
10000 CONTINUE
C
      RETURN
C
      ENTRY IF4COP
C
      NREC11=NOUTAU
      OPEN (UNIT=11,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC11)
      IRC=3
      WRITE (11,REC=IRC) RCCPL
      IF (TRCINT.OR.TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      ALLOCATE (IHELP(NOUTAU))
      IHELP=0
      JC=0
      DO K=1,NPTRGT
        DO J=1,10*NSTEP
          JC=JC+1
          IHELP(JC)=ICCPL1(J,K)
          IF (JC == NOUTAU) THEN
            IRC=IRC+1
            WRITE (11,REC=IRC) IHELP
            IF (TRCINT.OR.TRCFLE) 
     .        WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
            JC=0
          END IF
        END DO
      END DO
      IF (JC > 0) THEN
        IRC=IRC+1
        WRITE (11,REC=IRC) IHELP
        IF (TRCINT.OR.TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      END IF
      DEALLOCATE (IHELP)
      IRC=IRC+1
      WRITE (11,REC=IRC) ICCPL2
      IRC=IRC+1
      WRITE (11,REC=IRC) LCCPL
      IF (TRCINT.OR.TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
C
      IF (LSHORT) LSTOP=LSTP
C
      IF (.NOT.LSTOP) RETURN
C
      IF (.NOT.(LBALAN)) GOTO 11000
C
C  BALANCES, SHOULD BE DONE ONLY AT THE END OF BRAAMS RUN
C  AT THE END OF AN EIRENE RUN THE BALANCES MAY BE OFF AT LEAST AT
C  THE BEGINNING OF THE CYCLING PROCEDURE, BECAUSE THE PLASMA STILL
C  HAS TO ADJUST TO THE NEW SOURCES
C
C
C
11000 CONTINUE
C
      DEALLOCATE (NODBRICK)
      RETURN
C
      END


      subroutine kreuzprod (a,b,vnorm,c)
      use precision
      implicit none
      real(dp), intent(in) :: a(3), b(3), vnorm
      real(dp), intent(out) :: c(3)

      c(1) = (a(2)*b(3) - a(3)*b(2)) / vnorm
      c(2) = (a(3)*b(1) - a(1)*b(3)) / vnorm
      c(3) = (a(1)*b(2) - a(2)*b(1)) / vnorm

      return
      end subroutine kreuzprod
