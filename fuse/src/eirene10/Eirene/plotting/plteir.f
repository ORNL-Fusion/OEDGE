cdr  30.4.04:  call plttly for spectra corrected.
cdr            first bin (no.0) and last bin (no. nsts+1) contain
cdr            the fluxes outside the range of spectra.
cpb  30.7.04:  deal with switched off tallies
cdr  10.6.05:  further modifications of plot for spectra (text,
c              total, plot vs. wavelength, plot 2 spectra into same frame
c    7.12.06:  in call to rstrt: one argument was wrong: sgms_cop--> sgms_bgk
!pb  18.12.06: general checking of XMCP removed to allow plots of 
!              input tallies even is no Monte Carlo particle has been followed
!    10.01.07: ENTRY PLTEIR_REINIT added for reinitialization of EIRENE
C
C
      SUBROUTINE PLTEIR (ISTRA)
C
C  ISTRA IS THE STRATUM NUMBER. ISTRA=0 STANDS FOR: SUM OVER STRATA
C  PLOT PLASMA TALLIES ONLY ONCE.
C
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CESTIM
      USE CCONA
      USE CGRPTL
      USE CLOGAU
      USE CPLMSK
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CTRCEI
      USE CGEOM
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE COMSOU
      USE CTEXT
      USE COUTAU
      USE CSPEI

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA

      REAL(DP), ALLOCATABLE :: VECTOR(:,:),VECSAV(:,:),VSDVI(:,:)
      REAL(DP), ALLOCATABLE :: XSPEC(:),YSPEC(:,:),VSPEC(:,:),
     .          WLSPEC(:),YSPECWL(:,:),VSPECWL(:,:)
      REAL(DP) :: XXP3D_DUM(1), YYP3D_DUM(1)
      REAL(DP) :: YMN2(NPLT), YMX2(NPLT), YMNLG2(NPLT), YMXLG2(NPLT)
      REAL(DP) :: XMI, XMA, TMIN, TMAX, XI, XE, DEL, OUTAUI,
     .            SPCAN, SPC00,WL00,DE, DW
      INTEGER :: IR1(NPLT), IR2(NPLT), IRS(NPLT)
      INTEGER :: IXXE, IXXI, IYYE, IYYI, K, IINDEX, ISPC, NSPS, INULL,
     .           NF, NFT, I, IA, N, IXSET2, ISPZ, IALG, N1SDVI, ISAVE,
     .           IFIRST, IALV, ITL, JTAL, IBLD, ICURV, IE, IXSET3, IS,
     .           IERR, ICINC, IYSET3, IX, I2M, J, IRAD, I0, I1, I2, IT,
     .           INDX
      LOGICAL :: LPLOT2(NPLT), LSDVI(NPLT), LPLTT2, LINLOG, L_SAME
      CHARACTER(24) :: TXUNIT(NPLT), TXSPEC(NPLT)
      CHARACTER(24) :: TXUNT1, TXSPC1
      CHARACTER(72) :: TXTALL(NPLT)
      CHARACTER(72) :: TXTLL1
      CHARACTER(72) :: HEAD,  HEAD0, HEAD1, HEAD2, HEAD3, HEAD4,
     .                 HEAD5, HEAD6, HEAD7, HEAD8, HEAD9, HEAD10, TXHEAD
C
C
      SAVE
      DATA IFIRST/0/
      IF (IFIRST.EQ.0) THEN
        ISAVE=ISTRA
        IFIRST=1
      ENDIF
C
      IF (TRCPLT)
     .    WRITE (iunout,*) 'PLTEIR CALLED, ISTRA, XMCP: ',
     .                      ISTRA,XMCP(ISTRA)
C
C  NULLPUNKT AUF DEM PAPIER

      X0PL=10.
      Y0PL=3.
C  ACHSENLAENGEN
      LENX=25.
      LENY=20.
C  ACHSENUNTERTEILUNG VORGEGEBEN?
C  NEIN!
      STPSZX=0.
      STPSZY=0.
      INTNRX=0
      INTNRY=0
C  ACHSE LOGARITHMISCH?
      LOGX=.FALSE.
C     LOGY VIA INPUT
C  LOG. ACHSE MIN
      MINLY=0
C  LOG. ACHSE MAX
C     MAXLY WERDEN BERECHNET IN ANPSGL
C  ZEICHNE NETZLINIEN EIN
      GRIDX=.TRUE.
      GRIDY=.TRUE.
C  MACHE GRADE GRENZEN, X-ACHSE (Y ACHSE, NUR WENN TALZMI=TALZMA=666.)
      FITX=.TRUE.
C  NEW FRAME FOR EACH PICTURE IN PLTTLY
      L_SAME=.FALSE.
C
C
      IF (IESTR.EQ.ISTRA) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IF (XMCP(ISTRA).GT.1.) THEN
        IESTR=ISTRA
        CALL RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL SYMET(ESTIMV,NTALV,NRAD,NR1ST,NP2ND,NT3RD,
     .               NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
        ENDIF
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
        IF (XMCP(ISTRA).GT.1.) THEN
        IESTR=ISTRA
        CALL RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL SYMET(ESTIMV,NTALV,NRAD,NR1ST,NP2ND,NT3RD,
     .               NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
        ENDIF
      ELSE
        WRITE (iunout,*) 'ERROR IN PLTEIR: DATA FOR STRATUM ISTRA= ',
     .                    ISTRA
        WRITE (iunout,*) 'ARE NOT AVAILABLE. PLOTS ABANDONNED'
        RETURN
      ENDIF
C
      IF (ISTRA.EQ.0)
     .HEAD='SUM OVER STRATA
     .          '
      IF (ISTRA.NE.0) THEN
      HEAD='STRATUM NO.
     .          '
      WRITE (HEAD(13:15),'(I3)') ISTRA
      ENDIF
C
      HEAD0='VOLUME AVERAGED BACKGROUND TALLY, INPUT
     .           '
      HEAD1='DEFAULT VOLUME AVERAGED TALLY, TRACKLENGTH ESTIMATED
     .           '
      HEAD2='ADDITIONAL VOLUME AVERAGED TALLY, TRACKLENGTH ESTIMATED
     .           '
      HEAD3='ADDITIONAL VOLUME AVERAGED TALLY, COLLISION ESTIMATED
     .           '
      HEAD4='VOLUME AVERAGED TALLY, SNAPSHOT ESTIMATED
     .           '
      HEAD5='VOLUME AVERAGED TALLY, FOR COUPLING TO PLASMA CODE
     .           '
      HEAD6='BGK TALLY
     .           '
      HEAD7='ALGEBRAIC FUNCTION OF VOLUME AVERAGED TALLIES
     .           '
      HEAD8='RELATIVE STANDARD DEVIATION
     .           '
      HEAD9='SPECTRUM (VS. ENERGY, EV)
     .           '
      HEAD10='SPECTRUM (VS. WAVELENGTH, NM)
     .            '
C
      IALG=0
C
C  .......................................
C
C   LOOP OVER NVOLPL
C  .......................................
C
      IF (NVOLPL > 0) THEN
        ALLOCATE (VECTOR(NRAD,NPLT))
        IF (ANY(PLTL2D(1:NVOLPL).AND.PLTL3D(1:NVOLPL)))
     .     ALLOCATE (VECSAV(NRAD,NPLT))
        IF (ANY(PLTLER(1:NVOLPL))) THEN
          ALLOCATE (VSDVI(NRAD,NPLT))
          N1SDVI = NRAD
        ELSE
          ALLOCATE (VSDVI(1,1))
          N1SDVI = 1
        END IF
      END IF

      DO 10000 IBLD=1,NVOLPL
C
        IF (PLTL2D(IBLD).OR.PLTL3D(IBLD)) THEN
C
          DO 110 ICURV=1,NSPTAL(IBLD)
            JTAL=NPTALI(IBLD,ICURV)
C  REDO ALGEBRAIC TALLY IN CASE NFILEN=2 OR NFILEN=7
            IF (JTAL.EQ.NTALR.AND.IALG.EQ.0.AND.
     .          (NFILEN.EQ.2.OR.NFILEN.EQ.7)) THEN
              CALL ALGTAL
              IALG=1
              DO 105 IALV=1,NALVI
                CALL INTTAL (ALGV,VOLTAL,IALV,NALV,NSBOX_TAL,
     .                       ALGVI(IALV,ISTRA),
     .                       NR1TAL,NP2TAL,NT3TAL,NBMLT)
105           CONTINUE
            ENDIF
            ITL=IABS(JTAL)
C  PLOT OUTPUT TALLIES ONLY FOR STRATA WITH TWO OR MORE HISTORIES
            IF (JTAL.GT.0.AND.XMCP(ISTRA).LE.1) GOTO 10000
C  PLOT INPUT TALLIES ONLY ONCE PER ITERATION
            IF (JTAL.LT.0.AND.ISAVE.NE.ISTRA) GOTO 10000
            TXHEAD=HEAD0
            IF (JTAL.GT.0)     TXHEAD=HEAD1
            IF (JTAL.EQ.NTALA) TXHEAD=HEAD2
            IF (JTAL.EQ.NTALC) TXHEAD=HEAD3
            IF (JTAL.EQ.NTALT) TXHEAD=HEAD4
            IF (JTAL.EQ.NTALM) TXHEAD=HEAD5
            IF (JTAL.EQ.NTALB) TXHEAD=HEAD6
            IF (JTAL.EQ.NTALR) TXHEAD=HEAD7
            IF (TRCPLT) THEN
              CALL LEER(1)
              WRITE (iunout,*) 'PLOT REQUESTED FOR TALLY NO. ',JTAL
            ENDIF
C
            LPLTT2=.FALSE.
C
C .............................................
C
C  PUT TALLY ONTO ARRAY: VECTOR
C .............................................
C
C
            LSDVI(ICURV)=.FALSE.
            LPLOT2(ICURV)=.FALSE.
            ISPZ=ISPTAL(IBLD,ICURV)
            IF (JTAL.LT.0.) THEN
              NF=NFRSTP(ITL)
              DO 111 I=1,NRAD
111             VECTOR(I,ICURV)=0.
              IF (ISPZ.EQ.0) THEN
                SELECT CASE (ITL)
                CASE (1)
                  VECTOR(1:NSBOX,ICURV) = TEIN(1:NSBOX)
                CASE (2)
                  VECTOR(1:NSBOX,ICURV) = SUM(TIIN(1:NF,1:NSBOX),1)
                CASE (3)
                  VECTOR(1:NSBOX,ICURV) = DEIN(1:NSBOX)
                CASE (4)
                  VECTOR(1:NSBOX,ICURV) = SUM(DIIN(1:NF,1:NSBOX),1)
                CASE (5)
                  VECTOR(1:NSBOX,ICURV) = SUM(VXIN(1:NF,1:NSBOX),1)
                CASE (6)
                  VECTOR(1:NSBOX,ICURV) = SUM(VYIN(1:NF,1:NSBOX),1)
                CASE (7)
                  VECTOR(1:NSBOX,ICURV) = SUM(VZIN(1:NF,1:NSBOX),1)
                CASE (8)
                  VECTOR(1:NSBOX,ICURV) = BXIN(1:NSBOX)
                CASE (9)
                  VECTOR(1:NSBOX,ICURV) = BYIN(1:NSBOX)
                CASE (10)
                  VECTOR(1:NSBOX,ICURV) = BZIN(1:NSBOX)
                CASE (11)
                  VECTOR(1:NSBOX,ICURV) = BFIN(1:NSBOX)
                CASE (12)
                  VECTOR(1:NSBOX,ICURV) = SUM(ADIN(1:NF,1:NSBOX),1)
                CASE (13)
                  VECTOR(1:NSBOX,ICURV) = SUM(EDRIFT(1:NF,1:NSBOX),1)
                CASE (14)
                  VECTOR(1:NSBOX,ICURV) = VOL(1:NSBOX)
                CASE (15)
                  VECTOR(1:NSBOX,ICURV) = SUM(WGHT(1:NF,1:NSBOX),1)
                CASE (16)
                  VECTOR(1:NSBOX,ICURV) = BXPERP(1:NSBOX)
                CASE (17)
                  VECTOR(1:NSBOX,ICURV) = BYPERP(1:NSBOX)
                CASE DEFAULT
                  WRITE (iunout,*) ' WRONG TALLY NUMBER IN PLTEIR',
     .                        ' JTAL = ',JTAL
                  WRITE (iunout,*) ' NO PLOT PERFORMED '
                  CALL LEER(1)
                  GOTO 10000
                END SELECT
              ELSEIF (ISPZ.GT.0.AND.ISPZ.LE.NF) THEN
                SELECT CASE (ITL)
                CASE (1)
                  VECTOR(1:NSBOX,ICURV) = TEIN(1:NSBOX)
                CASE (2)
                  VECTOR(1:NSBOX,ICURV) = TIIN(MPLSTI(ISPZ),1:NSBOX)
                CASE (3)
                  VECTOR(1:NSBOX,ICURV) = DEIN(1:NSBOX)
                CASE (4)
                  VECTOR(1:NSBOX,ICURV) = DIIN(ISPZ,1:NSBOX)
                CASE (5)
                  VECTOR(1:NSBOX,ICURV) = VXIN(MPLSV(ISPZ),1:NSBOX)
                CASE (6)
                  VECTOR(1:NSBOX,ICURV) = VYIN(MPLSV(ISPZ),1:NSBOX)
                CASE (7)
                  VECTOR(1:NSBOX,ICURV) = VZIN(MPLSV(ISPZ),1:NSBOX)
                CASE (8)
                  VECTOR(1:NSBOX,ICURV) = BXIN(1:NSBOX)
                CASE (9)
                  VECTOR(1:NSBOX,ICURV) = BYIN(1:NSBOX)
                CASE (10)
                  VECTOR(1:NSBOX,ICURV) = BZIN(1:NSBOX)
                CASE (11)
                  VECTOR(1:NSBOX,ICURV) = BFIN(1:NSBOX)
                CASE (12)
                  VECTOR(1:NSBOX,ICURV) = ADIN(ISPZ,1:NSBOX)
                CASE (13)
                  VECTOR(1:NSBOX,ICURV) = EDRIFT(ISPZ,1:NSBOX)
                CASE (14)
                  VECTOR(1:NSBOX,ICURV) = VOL(1:NSBOX)
                CASE (15)
                  VECTOR(1:NSBOX,ICURV) = WGHT(ISPZ,1:NSBOX)
                CASE (16)
                  VECTOR(1:NSBOX,ICURV) = BXPERP(1:NSBOX)
                CASE (17)
                  VECTOR(1:NSBOX,ICURV) = BYPERP(1:NSBOX)
                CASE DEFAULT
                  WRITE (iunout,*) ' WRONG TALLY NUMBER IN PLTEIR',
     .                        ' JTAL = ',JTAL
                  WRITE (iunout,*) ' NO PLOT PERFORMED '
                  CALL LEER(1)
                  GOTO 10000
                END SELECT
              ELSE
                IF (TRCPLT) THEN
                  WRITE (iunout,*) 'SPECIES INDEX OUT OF RANGE '
                  WRITE (iunout,*) 'ICURV,ISPTAL(IBLD,ICURV) ',
     .                              ICURV,ISPZ
                  WRITE (iunout,*) 
     .              'ALL PLOTS FOR THIS TALLY TURNED OFF '
                ENDIF
                PLTL2D(IBLD)=.FALSE.
                PLTL3D(IBLD)=.FALSE.
                GOTO 110
              ENDIF
            ELSEIF (JTAL.GE.0) THEN
              IF (.NOT.LIVTALV(JTAL)) THEN
                WRITE (iunout,*) TXTTAL(1,JTAL)
                WRITE (iunout,*) 'TALLY SWITCHED OFF '
                WRITE (iunout,*) 'ALL PLOTS FOR THIS TALLY TURNED OFF '
                GOTO 10000
              END IF
              NFT=NFSTVI(ITL)
              NF=NFIRST(ITL)
              IF (ISPZ.EQ.0) THEN
                DO 121 I=1,NRAD
121               VECTOR(I,ICURV)=0.
                DO 122 K=1,NFT
                  DO 122 I=1,NRAD
                    VECTOR(I,ICURV)=VECTOR(I,ICURV)+
     .                              ESTIMV(NADDV(ITL)+K,NCLTAL(I))
122             CONTINUE
              ELSEIF (ISPZ.GT.0.AND.ISPZ.LE.NFT) THEN
                DO 125 I=1,NRAD
                  VECTOR(I,ICURV)=ESTIMV(NADDV(ITL)+ISPZ,NCLTAL(I))
125             CONTINUE
              ELSE
                IF (TRCPLT) THEN
                  WRITE (iunout,*) 'SPECIES INDEX OUT OF RANGE '
                  WRITE (iunout,*) 'ICURV,ISPTAL(IBLD,ICURV) ',
     .                              ICURV,ISPZ
                  WRITE (iunout,*) 
     .              'ALL PLOTS FOR THIS TALLY TURNED OFF '
                ENDIF
                PLTL2D(IBLD)=.FALSE.
                PLTL3D(IBLD)=.FALSE.
                GOTO 110
              ENDIF
C
              IF (PLTLER(IBLD)) THEN
C  CHECK IF STANDARD DEVIATION IS AVAILABLE FOR THIS TALLY
                DO 126 N=1,NSIGVI
                  IF (IIH(N).NE.JTAL) GOTO 126
                  IF (IGH(N).NE.ISPZ.AND.IGH(N).NE.0) GOTO 126
                  LSDVI(ICURV)=.TRUE.
                  DO 127 I=1,NRAD
                    VSDVI(I,ICURV)=SIGMA(N,NCLTAL(I))
127               CONTINUE
126             CONTINUE
              ENDIF
C
            ENDIF
C
            IF (PLTL2D(IBLD) .AND. PLTL3D(IBLD)) THEN
              DO 129 I=1,NRAD
                VECSAV(I,ICURV)=VECTOR(I,ICURV)
129           CONTINUE
            END IF
C
110       CONTINUE
C
C ...................................
C                                   .
C    VECTOR(IC,ICURV) IS SET NOW    .
C ...................................
C
C
          LOGY=PLTLLG(IBLD)
C
          IF (PLTL2D(IBLD)) THEN
C
C  ............................
C
C   SET ABSCISSA FOR 2D PLOT
C  ............................
C
C  SET ABSCISSA FROM GRID DATA, OR USE ONE OF THE INPUT OPTIONS:
            IXSET2=0
            IF (TALXMI(IBLD).NE.0..OR.TALXMA(IBLD).NE.0.) THEN
              XMI=TALXMI(IBLD)
              XMA=TALXMA(IBLD)
              IA= 100000000
              IE=-100000000
              DO ICURV=1,NSPTAL(IBLD)
                IA=MIN(IA,NPLIN2(IBLD,ICURV))
                IE=MAX(IE,NPLOT2(IBLD,ICURV))
              ENDDO
              DEL=IE-IA
              IF (XMI.LT.XMA) THEN
C  EQUIDISTANT IN LIN SCALE
                DO I=IA,IE
                  XXP2D(I)=XMI+(I-IA)/DEL*(XMA-XMI)
                ENDDO
              ELSEIF (XMI.GT.XMA.AND.XMI.GT.0.AND.XMA.GT.0) THEN
C  EQUIDISTANT IN LOG SCALE
                XI=LOG(XMA)
                XE=LOG(XMI)
                DO I=IA,IE
                  XXP2D(I)=EXP(XI+(I-IA)/DEL*(XE-XI))
                ENDDO
C  USER DEFINED ABSCISSA, XXP2D_USR
              ELSEIF (XMI.GT.XMA.AND.(XMI.LE.0.OR.XMA.LE.0)) THEN
                DO I=IA,IE
                  XXP2D(I)=XXP2D_USR(I,IBLD)
                ENDDO
              ENDIF
              IXSET2=1
              XMI=XXP2D(IA)*(1.+1.E-6)
              XMA=XXP2D(IE)/(1.+1.E-6)
              GOTO 139
            ENDIF
C
C  TRY DEFAULT OPTION TO SET PLOT GRID FROM 1.ST (RADIAL) GRID
C
            IXSET2=0
            IF (LEVGEO.EQ.1.OR.LEVGEO.EQ.2) THEN
C   USE RADIAL SURFACE CENTERED GRID "RHOSRF" FOR EACH POLOIDAL POSITION
              DO 130 I=1,NR1ST
                XXP2D(I)=RHOSRF(I)
130           CONTINUE
              DO 131 J=2,NP2ND*NT3RD*NBMLT
                DO 131 I=1,NR1ST
                  XXP2D(I+(J-1)*NR1ST)=XXP2D(I)
131            CONTINUE
              DO 138 I=NSURF+1,NRAD
138             XXP2D(I)=0.
              IXSET2=1
            ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
C   USE PERPEND. ARCLENGTH "BGLP" IN CASE OF POLYGON GRID AND NLPOL
              DO 133 I=1,NR1ST
                DO 133 J=1,NP2ND
                  DO 133 K=1,NT3RD
                    IRAD=I+((J-1)+(K-1)*NP2T3)*NR1P2
                    XXP2D(IRAD)=BGLP(I,J)
133           CONTINUE
              DO 136 I=NSURF+1,NRAD
136             XXP2D(I)=0.
              IXSET2=1
            ELSE
C   NO 2D PLOTOPTIONS AVAILABLE
            ENDIF
            XMI=XXP2D(NPLIN2(IBLD,1))*(1.+1.E-6)
            XMA=XXP2D(NPLOT2(IBLD,1))/(1.+1.E-6)
C
139         CONTINUE
C
            IF (IXSET2.NE.1) THEN
              WRITE (iunout,*) ' NO GRID SET FOR 2D PLOTTING '
              WRITE (iunout,*) ' IXSET2 = ',IXSET2
              WRITE (iunout,*) ' NO 2D PLOTTING IS DONE '
              GOTO 1000
            ENDIF
C
C  IN CASE OF LSMOT2, SET ZONE CENTERED ABSCISSA
C  GRID FROM SURFACE CENTERED  GRID  "X"
C
            IF (LSMOT2(IBLD)) THEN
              DO 137 J=1,NRAD-1
                XXP2D(J)=(XXP2D(J)+XXP2D(J+1))*0.5
137           CONTINUE
            ENDIF
C
            DO 140 ICURV=1,NSPTAL(IBLD)
              JTAL=NPTALI(IBLD,ICURV)
              ITL=IABS(JTAL)
              ISPZ=ISPTAL(IBLD,ICURV)
              IR1(ICURV)=NPLIN2(IBLD,ICURV)
              IR2(ICURV)=NPLOT2(IBLD,ICURV)
              IRS(ICURV)=NPLDL2(IBLD,ICURV)
C
              YMNLG2(ICURV)=1.D60
              YMXLG2(ICURV)=-1.D60
              IF (JTAL.GT.0.) THEN
                I0=0
                IF (NFRSTI(ITL).GT.1) I0=1
                INDX=NADDI(ITL)*NSTRAP+NFRSTI(ITL)*ISTRA+ISPZ+I0
                CALL FETCH_OUTAU (OUTAUI,JTAL,ISPZ,ISTRA,IUNOUT)
                IF (OUTAUI.EQ.0.) THEN
                  IF (TRCPLT) THEN
                    WRITE (iunout,*) 'TALLY NO. ',JTAL,
     .                               ' CURVE NO. ',ICURV
                    WRITE (iunout,*) 'NOT PLOTTED BECAUSE'
                    WRITE (iunout,*) 'ZERO INTEGRAL (OUTAU(INDX)=0.) '
                    WRITE (iunout,*) 'INDX,NADDI(JTAL),NFRSTI(JTAL),I0'
                    WRITE (iunout,*)  INDX,NADDI(ITL),NFRSTI(ITL),I0
                  ENDIF
                  YMN2(ICURV)=0.
                  YMX2(ICURV)=0.
                  YMNLG2(ICURV)=0.
                  YMXLG2(ICURV)=0.
                  GOTO 140
                ENDIF
              ENDIF
              LPLOT2(ICURV)=.TRUE.
              LPLTT2=.TRUE.
              I1=IR1(ICURV)
              I2=IR2(ICURV)
              I2M=I2-1
              IS=IRS(ICURV)
C
C YMNLG2, YMXLG2: REAL MAX/MIN, FOR LEGENDE ON 2D PLOT ONLY
              DO 141 I=I1,I2M,IS
                YMNLG2(ICURV)=MIN(YMNLG2(ICURV),VECTOR(I,ICURV))
141           CONTINUE
              DO 142 I=I1,I2M,IS
                YMXLG2(ICURV)=MAX(YMXLG2(ICURV),VECTOR(I,ICURV))
142           CONTINUE
C
C YMN2, YMX2: FOR AXIS
              FITY=.TRUE.
              IF (TALZMI(IBLD).NE.666.) THEN
                IF (.NOT.LOGY) FITY=.FALSE.
                YMN2(ICURV)=TALZMI(IBLD)
                DO 143 I=1,NRAD
                  VECTOR(I,ICURV)=MAX(YMN2(ICURV),VECTOR(I,ICURV))
143             CONTINUE
                IF (LOGY) YMN2(ICURV)=YMN2(ICURV)*(1.+1.E-6)
              ELSE
                YMN2(ICURV)=YMNLG2(ICURV)
              ENDIF
C
              IF (TALZMA(IBLD).NE.666.) THEN
                IF (.NOT.LOGY) FITY=.FALSE.
                YMX2(ICURV)=TALZMA(IBLD)
                DO 144 I=1,NRAD
                  VECTOR(I,ICURV)=MIN(YMX2(ICURV),VECTOR(I,ICURV))
144             CONTINUE
                IF (LOGY) YMX2(ICURV)=YMX2(ICURV)/(1.+1.E-6)
              ELSE
                YMX2(ICURV)=YMXLG2(ICURV)
              ENDIF
140         CONTINUE
C
C  PLOT ALL CURVES REQUESTED FROM THIS TALLY INTO ONE PICTURE
            DO 150 ICURV=1,NSPTAL(IBLD)
              JTAL=NPTALI(IBLD,ICURV)
              ITL=IABS(JTAL)
              ISPZ=ISPTAL(IBLD,ICURV)
              IF (ISPZ.EQ.0) THEN
                TXSPEC(ICURV)='SUM OVER SPECIES        '
                IF (JTAL.LT.0) TXUNIT(ICURV)=TXTPUN(1,ITL)
                IF (JTAL.LT.0) TXTALL(ICURV)=TXTPLS(1,ITL)
                IF (JTAL.GE.0) TXUNIT(ICURV)=TXTUNT(1,ITL)
                IF (JTAL.GE.0) TXTALL(ICURV)=TXTTAL(1,ITL)
              ELSE
                IF (JTAL.LT.0) THEN
                  TXTALL(ICURV)=TXTPLS(ISPZ,ITL)
                  TXSPEC(ICURV)=TXTPSP(ISPZ,ITL)
                  TXUNIT(ICURV)=TXTPUN(ISPZ,ITL)
                ELSE
                  TXTALL(ICURV)=TXTTAL(ISPZ,ITL)
                  TXSPEC(ICURV)=TXTSPC(ISPZ,ITL)
                  TXUNIT(ICURV)=TXTUNT(ISPZ,ITL)
                ENDIF
              ENDIF
150         CONTINUE
            IERR=0
            L_SAME=.FALSE.
            CALL PLTTLY (XXP2D,VECTOR,VSDVI,YMN2,YMX2,
     .             IR1,IR2,IRS,
     .             NSPTAL(IBLD),TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .             LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,LHIST2(IBLD),IERR,
     .             N1SDVI,NRAD,L_SAME)
            IF (TRCPLT) THEN
              IF (IERR.GT.0) THEN
                WRITE (iunout,*) '2D PLOT FOR TALLY NO. ',JTAL,
     .                           ' ABANDONED'
                WRITE (iunout,*) 'ERROR CODE FROM SUBR. PLTTLY: ',IERR
                WRITE (iunout,*) 'XMI,XMA ',XMI,XMA
                GOTO 1000
              ENDIF
              WRITE (iunout,*) '2D PLOT FOR TALLY NO. ',JTAL,' DONE'
              WRITE (iunout,*) 'XMIN= ',XMI,' XMAX= ',XMA
              DO 160 ICURV=1,NSPTAL(IBLD)
                IF (LPLOT2(ICURV))
     .            WRITE (iunout,*) 'ICURV= ',ICURV,
     .                        ' YMIN= ',YMNLG2(ICURV),
     .                        ' YMAX= ',YMXLG2(ICURV),
     .                        ' LSDVI= ',LSDVI(ICURV)
160           CONTINUE
            ENDIF
C
          ENDIF
C
1000      CONTINUE
C
C   3D PLOT GRID
C
          IF (PLTL3D(IBLD)) THEN
C
            DO 1040 ICURV=1,NSPTAL(IBLD)
              IF (PLTL2D(IBLD)) THEN
                DO 1035 I=1,NRAD
                   VECTOR(I,ICURV)=VECSAV(I,ICURV)
1035            CONTINUE
              END IF
C  SYMMETRY CONDITION AT POLAR ANGLE THETA=YIA AND THETA=2*PI+YIA
C  NOT READY: IXTL3 NOT DEFINED HERE. ENFORCE SYMMETRY AUTOMATICALLY EARLIER
C             IF (LEVGEO.EQ.2.AND.IYTL3.EQ.NP2ND) THEN
C               DO 1036 I=1,IXTL3
C1036             VECTOR(I+NP2NDM*NR1ST,ICURV)=VECTOR(I,ICURV)
C             ENDIF
1040        CONTINUE
C
C SET QUASIRECTANGULAR PLOT GRIDS XXP3D (IX), IX=1,IXTL3
C                             AND YYP3D (IY), IY=1,IYTL3
C FOR EACH 3D PICTURE
C
            IXSET3=0
            IYSET3=0
            IF (LEVGEO.EQ.1) THEN
              IF (NLRAD.AND..NOT.LPRAD3(IBLD)) THEN
                IXTL3=NR1ST
                DO 218 I=1,IXTL3
218               XXP3D(I)=RHOSRF(I)
                IXSET3=1
              ENDIF
              IF (NLTOR.AND.NLTRZ.AND..NOT.LPTOR3(IBLD)) THEN
                IYTL3=NT3RD
                DO 220 I=1,IYTL3
220               YYP3D(I)=ZSURF(I)
                DO I=1,NR1ST
                  DO J=1,NT3RD
                    XPOL(I,J)=RHOSRF(I)
                    YPOL(I,J)=ZSURF(J)
                  ENDDO
                ENDDO
                IYSET3=1
              ENDIF
              IF (NLPOL.AND..NOT.LPPOL3(IBLD)) THEN
                IYTL3=NP2ND
                DO 221 I=1,IYTL3
221               YYP3D(I)=PSURF(I)
                DO I=1,NR1ST
                  DO J=1,NP2ND
                    XPOL(I,J)=RHOSRF(I)
                    YPOL(I,J)=PSURF(J)
                  ENDDO
                ENDDO
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.2) THEN
C
              IF (NLRAD.AND..NOT.LPRAD3(IBLD)) THEN
                IXTL3=NR1ST
                DO 223 I=1,IXTL3
223               XXP3D(I)=RHOSRF(I)
                IXSET3=1
              ENDIF
              IF (NLTOR.AND.NLTRZ.AND..NOT.LPTOR3(IBLD)) THEN
                IYTL3=NT3RD
                DO 226 I=1,IYTL3
226               YYP3D(I)=ZSURF(I)
                IYSET3=1
              ENDIF
              IF (NLPOL.AND..NOT.LPPOL3(IBLD)) THEN
                IYTL3=NP2ND
                DO 225 I=1,IYTL3-1
225               YYP3D(I)=0.5*(PSURF(I+1)+PSURF(I))
                YYP3D(NP2ND)=PSURF(1)+PI2A
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.3) THEN
C
              IF (LPTOR3(IBLD)) THEN
                IXTL3=NR1ST
                DO 228 IX=1,IXTL3
228               XXP3D(IX)=IX
                IXSET3=1
                IYTL3=NP2ND
                DO 230 IX=1,IYTL3
230               YYP3D(IX)=IX
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.4) THEN
C
              IF (LPTOR3(IBLD)) THEN
                IXSET3=1
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.5) THEN
C
!              IF (LPTOR3(IBLD)) THEN
                IXSET3=1
                IYSET3=1
!              ENDIF
            ELSE
              WRITE (iunout,*) '3D PLOT OPTION TO BE WRITTEN, LEVGEO '
              CALL EXIT_OWN(1)
            ENDIF
C
C  LOOP ICURV=1,....
C
            ICINC=1
            IF (LVECT3(IBLD).OR.LRPVC3(IBLD)) ICINC=2
            DO 1160 ICURV=1,NSPTAL(IBLD),ICINC
              JTAL=NPTALI(IBLD,ICURV)
              ITL=IABS(JTAL)
              ISPZ=ISPTAL(IBLD,ICURV)
              IF (ISPZ.EQ.0) THEN
                TXSPEC(1)='SUM OVER SPECIES        '
                IF (JTAL.LT.0) TXUNT1=TXTPUN(1,ITL)
                IF (JTAL.LT.0) TXTLL1=TXTPLS(1,ITL)
                IF (JTAL.GE.0) TXUNT1=TXTUNT(1,ITL)
                IF (JTAL.GE.0) TXTLL1=TXTTAL(1,ITL)
              ELSE
                IF (JTAL.LT.0) THEN
                  TXTLL1=TXTPLS(ISPZ,ITL)
                  TXSPC1=TXTPSP(ISPZ,ITL)
                  TXUNT1=TXTPUN(ISPZ,ITL)
                ELSE
                  TXTLL1=TXTTAL(ISPZ,ITL)
                  TXSPC1=TXTSPC(ISPZ,ITL)
                  TXUNT1=TXTUNT(ISPZ,ITL)
                ENDIF
              ENDIF
C
              IF (IXSET3+IYSET3.LT.2) THEN
                WRITE (iunout,*) ' NO GRIDS SET FOR 3D PLOTTING '
                WRITE (iunout,*) ' IXSET3,IYSET3 = ',IXSET3,IYSET3
                WRITE (iunout,*) ' NO 3D PLOTTING IS DONE '
                GOTO 10000
              ENDIF
C
              LINLOG=PLTLLG(IBLD)
              TMIN=TALZMI(IBLD)
              TMAX=TALZMA(IBLD)
C
1200          CONTINUE
C
C
C  CONTOUR PLOTS
              IF (LCNTR3(IBLD)) THEN
                CALL ISOLNE (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  WRITE FILES FOR RAPS PLOTS
              ELSEIF (LRAPS3(IBLD)) THEN
                CALL RPSCOL (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D_DUM,YYP3D_DUM,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  3D HISTOGRAM
              ELSEIF (LHIST3(IBLD)) THEN
                CALL PL3DPG (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,TALW1(IBLD),TALW2(IBLD),
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  3D SURFACE PLOTS, IN CUBE
              ELSEIF (LSMOT3(IBLD)) THEN
                CALL PLOT3D (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,TALW1(IBLD),TALW2(IBLD),
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  VECTOR FIELD PLOT
              ELSEIF (LVECT3(IBLD)) THEN
                IXXI=NPLI13(IBLD,ICURV)
                IXXE=NPLO13(IBLD,ICURV)
                IYYI=NPLI23(IBLD,ICURV)
                IYYE=NPLO23(IBLD,ICURV)
                CALL VECLNE (VECTOR(1,ICURV),
     .                       VECTOR(1,ICURV+1),IBLD,ICURV,
     .                       IXXI,IXXE,IYYI,IYYE,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  WRITE FILES FOR RAPS VECTOR PLOTS
              ELSEIF (LRPVC3(IBLD)) THEN
                CALL RPSVEC (VECTOR(1,ICURV),
     .                       VECTOR(1,ICURV+1),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D_DUM,YYP3D_DUM,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
              ELSE
                WRITE (iunout,*) 'NO 3D PLOT OPTION FOR IBLD= ',IBLD
                GOTO 10000
              ENDIF
C
C   PLOT STANDARD DEVIATION PROFILE FOR THIS TALLY, IF REQUESTED
C
              IF (LVECT3(IBLD).OR.LRPVC3(IBLD)) GOTO 1160
              IF (PLTLER(IBLD).AND.LSDVI(ICURV)) THEN
                TXUNT1='%                       '
                TXHEAD=HEAD8
                INULL=0
                DO 1222 I=1,NRAD
                  VECTOR(I,ICURV)=VSDVI(I,ICURV)
                  IF (LRAPS3(IBLD).AND.
     .               (ABS(VECTOR(I,ICURV)) < EPS30)) THEN
                    VECTOR(I,ICURV)=101._DP
                    INULL = INULL + 1
                  END IF
1222            CONTINUE
                LINLOG=.FALSE.
                TMIN=0.
                TMAX=100.
                IF (LRAPS3(IBLD).AND.INULL.GT.0) THEN
                  TMAX=101.
                  WRITE (IUNOUT,*) 'RAPS GRAPHICS FOR STD. DEVIATION:  '
                  WRITE (IUNOUT,*) 'IBLD, ICURV ',IBLD, ICURV
                  WRITE (IUNOUT,*)  INULL, ' CELLS WITH 0 HISTORIES '
                  WRITE (IUNOUT,*) 'STD. DEV. SET = 101% IN THESE CELLS'
                  WRITE (IUNOUT,*) 'TO PERMIT SPECIAL CHOICE OF COLOUR '
                END IF
                LSDVI(ICURV)=.FALSE.
                GOTO 1200
              ENDIF
C
1160        CONTINUE
C  LOOP ICURV FINISHED
          ENDIF
C
        ELSEIF (NSPTAL(IBLD).GT.0) THEN
          WRITE (iunout,*) 'PLOT REQUEST FOR TALLY NO. ',JTAL,
     .                     ' BUT NEITHER'
          WRITE (iunout,*) 
     .       'PLTL2D NOR PLTL3D TRUE. NO PLOT FOR THIS TALLY'
        ENDIF
C
10000 CONTINUE
C  LOOP IBLD FINISHED
C

      DO ISPC=1,NADSPC
C  THERE ARE NSPS BINS, AND NSPS+1 ENERGY BIN BOUNDARIES
        NSPS=ESTIML(ISPC)%PSPC%NSPC
        ALLOCATE (XSPEC(NSPS+1))
        ALLOCATE (YSPEC(NSPS+1,1))
        ALLOCATE (VSPEC(NSPS+1,1))
        SPCAN=ESTIML(ISPC)%PSPC%SPCMIN
        SPC00=ESTIML(ISPC)%PSPC%ESP_00
C  x axis: cell faces
        DO I=1,NSPS+1
          XSPEC(I)=SPCAN+(I-1)*ESTIML(ISPC)%PSPC%SPCDEL-SPC00
        END DO
C  y axis: cell averages (approx: cell centres)
        DO I=1,NSPS
          YSPEC(I,1)=ESTIML(ISPC)%PSPC%SPC(I)
          IF (NSIGI_SPC > 0) VSPEC(I,1)=ESTIML(ISPC)%PSPC%SGM(I)
        END DO

        YMN2(1)=MINVAL(YSPEC(1:NSPS,1))
        YMX2(1)=MAXVAL(YSPEC(1:NSPS,1))
        IF (ABS(YMX2(1)-YMN2(1)) < EPS30) YMX2(1) = YMN2(1) + 1._dp
        YMNLG2(1)=YMN2(1)
        YMXLG2(1)=YMX2(1)
        LSDVI(1)=NSIGI_SPC > 0
        LPLOT2(1)=.TRUE.
        IR1(1)=1
        IR2(1)=NSPS+1
        IRS(1)=1
        XMI=XSPEC(1)
        XMA=XSPEC(NSPS+1)
        LOGY=.FALSE.
        FITY=.FALSE.
        IF (ESTIML(ISPC)%PSPC%ISRFCLL == 0) THEN
         TXTALL(1)='SPECTRUM FOR SURFACE        PARTICLE TYPE        '//
     .            'SPECIES                '
        ELSE
         TXTALL(1)='SPECTRUM FOR CELL           PARTICLE TYPE        '//
     .            'SPECIES                '
        ENDIF
        WRITE (TXTALL(1)(22:27),'(I6)') ESTIML(ISPC)%PSPC%ISPCSRF
        WRITE (TXTALL(1)(43:48),'(I6)') ESTIML(ISPC)%PSPC%IPRTYP
        WRITE (TXTALL(1)(58:63),'(I6)') ESTIML(ISPC)%PSPC%IPRSP
        IT = ESTIML(ISPC)%PSPC%ISPCTYP
        TXSPEC=REPEAT(' ',24)
        TXUNIT=REPEAT(' ',24)
        IF (IT == 1) TXUNIT='AMP/BIN(EV)             '
        IF (IT == 2) TXUNIT='WATT/BIN(EV)            '
        TXHEAD=REPEAT(' ',72)
        TXHEAD(1:30)=HEAD9(1:30)
        TXHEAD(32:42)='INTEGRAL: '
        WRITE (TXHEAD(43:55),'(ES12.4)') ESTIML(ISPC)%PSPC%SPCINT
        IERR=0
        L_SAME=.TRUE.
        IF (ISPC.EQ.1) L_SAME=.FALSE.
        L_SAME=ESTIML(ISPC)%PSPC%SPC_SAME .NE. 1.D0
        CALL PLTTLY (XSPEC,YSPEC,VSPEC,YMN2,YMX2,
     .       IR1,IR2,IRS,
     .       1,TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .       LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,.TRUE.,IERR,
     .       NSPS+1,NSPS+1,L_SAME)
        DEALLOCATE (XSPEC)
        DEALLOCATE (YSPEC)
        DEALLOCATE (VSPEC)

      END DO

C  NOW REPEAT SAME PLOTS, BUT VS. WAVELENGTH

      DO ISPC=1,NADSPC
C  THERE ARE NSPS BINS, AND NSPS+1 ENERGY BIN BOUNDARIES
        NSPS=ESTIML(ISPC)%PSPC%NSPC
        ALLOCATE (XSPEC(NSPS+1))
        ALLOCATE (YSPEC(NSPS+1,1))
        ALLOCATE (VSPEC(NSPS+1,1))
        IF (NPHOTI > 0) THEN
          ALLOCATE (WLSPEC(NSPS+1))
          ALLOCATE (YSPECWL(NSPS+1,1))
          ALLOCATE (VSPECWL(NSPS+1,1))
        END IF
        SPCAN=ESTIML(ISPC)%PSPC%SPCMIN
        SPC00=ESTIML(ISPC)%PSPC%ESP_00
C  x axis: cell faces
        DO I=1,NSPS+1
          XSPEC(I)=SPCAN+(I-1)*ESTIML(ISPC)%PSPC%SPCDEL
        END DO
C  y axis: cell averages (approx: cell centres)
        DO I=1,NSPS
          YSPEC(I,1)=ESTIML(ISPC)%PSPC%SPC(I)
          IF (NSIGI_SPC > 0) VSPEC(I,1)=ESTIML(ISPC)%PSPC%SGM(I)
        END DO

        IF (NPHOTI > 0) THEN
C  PLOT ALSO VS. WAVELENGTH (NM)
          WL00            =HPCL/MAX(1.E-6_DP,SPC00)*1.E7_DP
C  x axis: cell faces
          DO I=1,NSPS+1
            WLSPEC(NSPS-I+1+1)=HPCL/MAX(1.E-6_DP,XSPEC(I))*1.E7_DP
            WLSPEC(NSPS-I+1+1)=WLSPEC(NSPS-I+1+1)-WL00
          END DO
C  y axis: cell averages (approx: cell centres)
          DO I=1,NSPS
            YSPECWL(NSPS-I+1,1)=YSPEC(I,1)
            IF (NSIGI_SPC > 0) VSPECWL(NSPS-I+1,1)=VSPEC(I,1)
          END DO
C  rescaling:  flux/ev to flux/nm
          DO I=1,NSPS
            DE=XSPEC(NSPS-I+1+1)-XSPEC(NSPS-I+1)
            DW=WLSPEC(I+1)-WLSPEC(I)
            YSPECWL(I,1) = YSPECWL(I,1)*DE/DW
          END DO
        END IF

        DEALLOCATE (XSPEC)
        DEALLOCATE (YSPEC)
        DEALLOCATE (VSPEC)

        IF (NPHOTI > 0) THEN
          YMN2(1)=MINVAL(YSPECWL(1:NSPS,1))
          YMX2(1)=MAXVAL(YSPECWL(1:NSPS,1))
          IF (ABS(YMX2(1)-YMN2(1)) < EPS30) YMX2(1) = YMN2(1) + 1._dp
          YMNLG2(1)=YMN2(1)
          YMXLG2(1)=YMX2(1)
          LSDVI(1)=NSIGI_SPC > 0
          LPLOT2(1)=.TRUE.
          IR1(1)=1
          IR2(1)=NSPS+1
          IRS(1)=1
          XMI=WLSPEC(1)
          XMA=WLSPEC(NSPS+1)
          LOGY=.TRUE.
          FITY=.TRUE.
          IF (ESTIML(ISPC)%PSPC%ISRFCLL == 0) THEN
            TXTALL(1)=
     .      'SPECTRUM FOR SURFACE        PARTICLE TYPE        '//
     .      'SPECIES                '
          ELSE
            TXTALL(1)=
     .      'SPECTRUM FOR CELL           PARTICLE TYPE        '//
     .      'SPECIES                '
          END IF
          WRITE (TXTALL(1)(22:27),'(I6)') ESTIML(ISPC)%PSPC%ISPCSRF
          WRITE (TXTALL(1)(43:48),'(I6)') ESTIML(ISPC)%PSPC%IPRTYP
          WRITE (TXTALL(1)(58:63),'(I6)') ESTIML(ISPC)%PSPC%IPRSP
          TXSPEC=REPEAT(' ',24)
          TXUNIT=REPEAT(' ',24)
          IF (IT == 1) TXUNIT='AMP/BIN(NM),            '
          IF (IT == 2) TXUNIT='WATT/BIN(NM),           '
          TXHEAD=REPEAT(' ',72)
          TXHEAD(1:30)=HEAD10(1:30)
          TXHEAD(32:42)='INTEGRAL: '
          WRITE (TXHEAD(43:55),'(ES12.4)') ESTIML(ISPC)%PSPC%SPCINT
          IERR=0
          L_SAME=.TRUE.
          IF (ISPC.EQ.1) L_SAME=.FALSE.
          L_SAME=ESTIML(ISPC)%PSPC%SPC_SAME .NE. 1.D0
          CALL PLTTLY (WLSPEC,YSPECWL,VSPECWL,YMN2,YMX2,
     .       IR1,IR2,IRS,
     .       1,TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .       LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,.TRUE.,IERR,
     .       NSPS+1,NSPS+1,L_SAME)
          DEALLOCATE (WLSPEC)
          DEALLOCATE (YSPECWL)
          DEALLOCATE (VSPECWL)
        END IF
      END DO


      IF (ALLOCATED(VECTOR)) DEALLOCATE(VECTOR)
      IF (ALLOCATED(VECSAV)) DEALLOCATE(VECSAV)
      IF (ALLOCATED(VSDVI))  DEALLOCATE(VSDVI)
      RETURN

C     the following ENTRY is for reinitialization of EIRENE (DMH)
      
      ENTRY PLTEIR_REINIT
      IFIRST = 0
      return
      END
