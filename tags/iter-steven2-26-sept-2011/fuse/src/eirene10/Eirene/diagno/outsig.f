C
C
C*DK OUTSIG
      SUBROUTINE OUTSIG(ENSAVE,L_CHOR)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMPRT, ONLY: IUNOUT
      USE CLOGAU
      USE CPLOT
      USE COMSIG
      USE CTEXT
      USE PHOTON
      USE CPLMSK

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ENSAVE(NCHOR,NCHEN)
      LOGICAL, INTENT(IN) :: L_CHOR(NCHOR)

      REAL(DP) :: XPLEN(NCHEN)
      REAL(DP) :: DUMFFD(NCHEN),WLDUMFFD(NCHEN),
     .            DELENE(NCHEN),DELWL(NCHEN),
     .            DEL_HELP(NCHEN), VSPEC(NCHEN)
      REAL(DP) :: YMN2(1), YMX2(1), YMNLG2(1), YMXLG2(1), XMI, XMA
      REAL(DP) :: TIMA, DUMTIL, AH, TIMI, DEL1, DEL2, SUM,WLSHFT
      INTEGER :: ICURV, ICHORI, I, NCHNI, ISK, JSK, IERR
      INTEGER :: IR1(1), IR2(1), IRS(1)
      LOGICAL :: LPLOT2(1), LSDVI(1), L_SAME
      CHARACTER(8) :: TSAFE, TEXTSS
      DATA TEXTSS /'SUM     '/
      CHARACTER(72) :: TXHEAD, TXTALL(1)
      CHARACTER(24) :: TXUNIT(1), TXSPEC(1)
C
      CALL PAGE
      NCHNI=IABS(NCHENI)
      TSAFE=TEXTS(1)
      ICURV=0
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

      DO 100  ICHORI=1,NCHORI
        IF (.NOT.L_CHOR(ICHORI)) GOTO 100
C
        CALL LEER(2)
        WRITE (iunout,*) 'NUMBER OF DETECTOR = ',ICHORI
        WRITE (iunout,*) TXTSIG(ICHORI)
        CALL LEER(1)
        TXHEAD=REPEAT(' ',72)
        IF (NCHTAL(ICHORI).EQ.1) THEN
          CALL HEADNG('CX-DETECTOR SIGNAL, MAXW. NEUTRL. DISTR. ',41)
          TXHEAD(1:41) = 
     .      'CX-DETECTOR SIGNALS, MAXW. NEUTRL. DISTR.'
        ELSEIF (NCHTAL(ICHORI).EQ.2) THEN
          CALL HEADNG('H ALPHA-DETECTOR SIGNAL: #/S/CM2/STERAD  ',41)
          TXHEAD(1:41) = 
     .      'H ALPHA-DETECTOR SIGNALS: #/S/CM2/STERAD '
        ELSEIF (NCHTAL(ICHORI).EQ.3) THEN
          CALL HEADNG('SIDE ON SPECTRAL RADIANCE: #/S/CM2/NM/STERAD',44)
          TXHEAD(1:44) =
     .      'SIDE ON SPECTRAL RADIANCE                   '
        ELSEIF (NCHTAL(ICHORI).EQ.4) THEN
          CALL HEADNG('DETECTOR SIGNALS VIA SPECTRA ALONG LOS',38)
          TXHEAD(1:44) =
     .      'DETECTOR SIGNALS VIA SPECTRA ALONG LOS'
        ELSEIF (NCHTAL(ICHORI).EQ.10) THEN
          CALL HEADNG('USER DEFINED LINE INTEGRAL (SUBR. SIGUSR)',41)
          TXHEAD(1:41) = 
     .      'USER DEFINED LINE INTEGRAL (SUBR. SIGUSR)'
        ELSE
          CALL HEADNG('LINE INTEGRAL TESTING OPTION: NO SIGNALS)',41)
          TXHEAD(1:41) = 
     .      'LINE INTEGRAL TESTING OPTION: NO SIGNALS)'
        ENDIF
        CALL LEER(2)
C
        TXTALL(1) = TXTSIG(ICHORI)
        CALL LEER(1)
        CALL MASR3 ('1ST POINT,"PIVOT POINT" ',
     .               XPIVOT(ICHORI),YPIVOT(ICHORI),ZPIVOT(ICHORI))
        CALL MASR3 ('2ND POINT, INSIDE VOLUME',
     .               XCHORD(ICHORI),YCHORD(ICHORI),ZCHORD(ICHORI))
        DO 30 I=1,NCHNI
          ENERGY(I)=ENSAVE(ICHORI,I)
30      CONTINUE
C...................................................
C  CX SPECTRA, ATOMS.  100 -- 199
C...................................................
        IF ((NCHTAL(ICHORI).EQ.1) .OR. (NCHTAL(ICHORI).EQ.4)) THEN
          IF (NSPSPZ(ICHORI).EQ.0) THEN
            TEXTS(1)=TEXTSS
          ENDIF

          IF (NCHTAL(ICHORI).EQ.4) THEN
            DO I=2,NCHNI-1
              DEL1=ENERGY(I+1)-ENERGY(I)
              DEL2=ENERGY(I)-ENERGY(I-1)
!pb              DELENE(I)=0.5_DP*(DEL1+DEL2)
              DELENE(I)=DEL1
            ENDDO
            IF (NCHNI > 1) THEN
C  FIRST INTERVAL (HALF SIZE)
!pb              DELENE(1)=0.5*(ENERGY(2)-ENERGY(1))
              DELENE(1)=ENERGY(2)-ENERGY(1)
C  LAST INTERVAL (HALF SIZE)
!pb              DELENE(NCHNI)=0.5*(ENERGY(NCHNI)-ENERGY(NCHNI-1))
              DELENE(NCHNI)=ENERGY(NCHNI)-ENERGY(NCHNI-1)

            ELSE
              DELENE(1) = 0._DP
            END IF
          END IF
C
          DUMTIL=TILINE(ICHORI)
          SUM=0._DP
          DO 70 I=1,NCHNI
            AH=FUFFER(ICHORI,I)
            DUMFFD(I)=MAX(1.E-10_DP,AH)
            IF (NCHTAL(ICHORI).EQ.4) SUM=SUM+AH*DELENE(I)
            XPLEN(I)=ENERGY(I)
70        CONTINUE

          CALL MASRR2('ENERGY,CXFLUX         ',ENERGY,DUMFFD,NCHNI)
          IF (NCHTAL(ICHORI).EQ.4) THEN
            CALL LEER(1)
            CALL MASR1('INTEGR. ',SUM)
            CALL LEER(1)
          END IF
          CALL MASR1('INP. TEM.',TINP(ICHORI))
          CALL MASR1('DT. TMP.',DUMTIL)
          CALL MASAGE ('FITTING RANGE:  TIMIN,TIMAX=                 ')
          TIMI=NSPINI(ICHORI)*TINP(ICHORI)
          TIMA=NSPEND(ICHORI)*TINP(ICHORI)
          CALL MASR2 ('TIMIN,TIMAX=    ',TIMI,TIMA)
          CALL LEER(2)
C
C  PREPARE DATA FOR PLOT OF SPECTRUM NO ICHORI
          IF (PLSPEC) THEN
            L_SAME = NSPNEW(ICHORI).NE.1
C  INITALIZE NEW PICTURE
            IF (NSPNEW(ICHORI).EQ.1) THEN
              IF (NSPSCL(ICHORI).EQ.0) THEN
                LOGX=.FALSE.
                LOGY=.FALSE.
              ELSEIF (NSPSCL(ICHORI).EQ.1) THEN
                LOGX=.FALSE.
                LOGY=.TRUE.
              ELSEIF (NSPSCL(ICHORI).EQ.2) THEN
                LOGX=.TRUE.
                LOGY=.FALSE.
              ELSEIF (NSPSCL(ICHORI).EQ.3) THEN
                LOGX=.TRUE.
                LOGY=.TRUE.
              ENDIF
            ENDIF
C  PLOT
            YMN2(1)=MINVAL(DUMFFD(1:NCHNI))
            YMX2(1)=MAXVAL(DUMFFD(1:NCHNI))
            IF (ABS(YMX2(1)-YMN2(1)) < EPS30) YMX2(1) = YMN2(1) + 1._dp
            YMNLG2(1)=YMN2(1)
            YMXLG2(1)=YMX2(1)
            LSDVI(1)=.FALSE.
            LPLOT2(1)=.TRUE.
            IR1(1)=1
            IR2(1)=NCHNI
            IRS(1)=1
            XMI=ENERGY(1)
            XMA=ENERGY(NCHNI)
            FITY=.FALSE.
            TXSPEC=REPEAT(' ',24)
            TXUNIT=REPEAT(' ',24)
            ICURV=ICURV+1
            CALL PLTTLY (ENERGY,DUMFFD,VSPEC,YMN2,YMX2,
     .           IR1,IR2,IRS,
     .           1,TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .           LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,.TRUE.,IERR,
     .           NCHNI,NCHNI,L_SAME)
          ENDIF
C
          TEXTS(1)=TSAFE
C..............................................
C  PHOTONS, LINE INTENSITY.  200-- 299
C..............................................
        ELSEIF (NCHTAL(ICHORI).EQ.2) THEN
C
          DUMTIL=FUFFER(ICHORI,1)
          CALL MASR1('H ALPHA ',DUMTIL)
          CALL LEER(2)

C...............................................
C  PHOTONS, SIDE ON SPECTRA.  300-- 399
C...............................................
        ELSEIF (NCHTAL(ICHORI).EQ.3) THEN

          IF (NSPSPZ(ICHORI).EQ.0) THEN
            TEXTS(1)=TEXTSS
          ENDIF

          DO I=2,NCHNI-1
            DEL1=ENERGY(I+1)-ENERGY(I)
            DEL2=ENERGY(I)-ENERGY(I-1)
            DELENE(I)=0.5_DP*(DEL1+DEL2)
          ENDDO
          IF (NCHNI > 1) THEN
C  FIRST INTERVAL (HALF SIZE)
            DELENE(1)=0.5*(ENERGY(2)-ENERGY(1))
C  LAST INTERVAL (HALF SIZE)
            DELENE(NCHNI)=0.5*(ENERGY(NCHNI)-ENERGY(NCHNI-1))

          ELSE
            DELENE(1) = 0._DP
          END IF

C  INTEGRATE SPECTRA, IN ENERGY PICTURE: SUM
C  SET ORDINATES, FOR PLOT AND PRINTOUT: DUMFFD
          SUM=0._DP
          DO 310 I=1,NCHNI
            AH=FUFFER(ICHORI,I)
            DUMFFD(NCHNI-I+1)=MAX(1.E-30_DP,AH)
            SUM=SUM+AH*DELENE(I)
c           write (iunout,*) 'outsig ',delene(i),energy(i),AH,Sum
310       CONTINUE

C  CONVERT FROM ENERGY TO WAVELENGTH

C  1.) INVERT SCALE, FOR CONVERSION TO WAVELENGTH:
          DEL_HELP=DELENE
          DO I=1,NCHNI
            DELENE(I)=DEL_HELP(NCHNI-I+1)
            DUMFFD(NCHNI-I+1)=MAX(1.E-30_DP,FUFFER(ICHORI,I))
          ENDDO
C  2.)  CONVERT FROM EV TO NM
          WLSHFT=0._DP
          IF (ABS(ESHIFT(ICHORI)).GT.1E-20)
     .            WLSHFT=HPCL/ESHIFT(ICHORI)*1.D7
          DO 370 I=1,NCHNI
            XPLEN(NCHNI-I+1)=HPCL/ENERGY(I)*1.D7-WLSHFT
370       CONTINUE

          DO I=2,NCHNI-1
            DEL1=XPLEN(I+1)-XPLEN(I)
            DEL2=XPLEN(I)-XPLEN(I-1)
            DELWL(I)=0.5_DP*(DEL1+DEL2)
          ENDDO
C  FIRST INTERVAL (HALF SIZE)
          DELWL(1)=0.5*(XPLEN(2)-XPLEN(1))
C  LAST INTERVAL (HALF SIZE)
          DELWL(NCHNI)=0.5*(XPLEN(NCHNI)-XPLEN(NCHNI-1))

C  3.)  CONVERT SPECTRAL DENSITY FROM EV TO NM
          WLDUMFFD=DUMFFD*DELWL/DELENE
C
          CALL MASRR2('WAVEL. ,RADIATIVE FLUX ',XPLEN,WLDUMFFD,NCHNI)
          CALL MASR1('INTEGR. ',SUM)
          CALL LEER(2)
C
C  PREPARE DATA FOR PLOT OF SPECTRUM NO ICHORI
          IF (PLSPEC) THEN
            L_SAME = NSPNEW(ICHORI).NE.1
C  INITALIZE NEW PICTURE
            IF (NSPNEW(ICHORI).EQ.1) THEN
              IF (NSPSCL(ICHORI).EQ.0) THEN
                LOGX=.FALSE.
                LOGY=.FALSE.
              ELSEIF (NSPSCL(ICHORI).EQ.1) THEN
                LOGX=.FALSE.
                LOGY=.TRUE.
              ELSEIF (NSPSCL(ICHORI).EQ.2) THEN
                LOGX=.TRUE.
                LOGY=.FALSE.
              ELSEIF (NSPSCL(ICHORI).EQ.3) THEN
                LOGX=.TRUE.
                LOGY=.TRUE.
              ENDIF
            ENDIF
C  PLOT, IN WAVELENGTH SCALE
            YMN2(1)=MINVAL(WLDUMFFD(1:NCHNI))
            YMX2(1)=MAXVAL(WLDUMFFD(1:NCHNI))
            IF (ABS(YMX2(1)-YMN2(1)) < EPS30) YMX2(1) = YMN2(1) + 1._dp
            YMNLG2(1)=YMN2(1)
            YMXLG2(1)=YMX2(1)
            LSDVI(1)=.FALSE.
            LPLOT2(1)=.TRUE.
            IR1(1)=1
            IR2(1)=NCHNI
            IRS(1)=1
            XMI=XPLEN(1)
            XMA=XPLEN(NCHNI)
            FITY=.FALSE.
            TXSPEC=REPEAT(' ',24)
            TXUNIT=REPEAT('PHOTONS/S/CM2/NM/STERAD ',24)
            ICURV=ICURV+1
            CALL PLTTLY (XPLEN,WLDUMFFD,VSPEC,YMN2,YMX2,
     .           IR1,IR2,IRS,
     .           1,TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .           LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,.FALSE.,IERR,
     .           NCHNI,NCHNI,L_SAME)
          ENDIF
C
          TEXTS(1)=TSAFE
C

        ELSEIF (NCHTAL(ICHORI).EQ.10) THEN
          write (iunout,*) 'printout for user defined line integral '
          write (iunout,*) 'still to be written in subr. outsig '
        ENDIF

100   CONTINUE
      CALL LEER(2)

      RETURN
      END
