c april05:  *sqrt(ze) moved from here (for cx spectra) into sigcx
c april06:  restriction to iphot.eq.isp in case of los-radiances
C
C
      SUBROUTINE EIRENE_SGNAL(ICHORI,IISTR,ISP,LCHOR)
C
C  THIS SUBROUTINE CALCULATES LINE INTEGRATED SIGNALS, USING THE EIRENE
C  VOLUME AVERAGED TALLIES AND THE PLASMA BACKGROUND DATA.
C  THERE MAY BE A CONTRIBUTION DIRECTLY FROM A PRIMARY SOURCE,
C  DUE TO DIRECT EMISSION FROM SOURCE INTO LINE OF SIGHT,
C  AS WELL AS A SECONDARY SOURCE (POST COLLISION) CONTRIBUTION, DUE TO
C  SCATTERING INTO THE LINE OF SIGHT
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPLOT
      USE EIRMOD_COMSIG
      USE EIRMOD_CGRID
      USE EIRMOD_CTRCEI
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CTEXT
      USE EIRMOD_CUPD
 
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: ICHORI, IISTR, ISP
      REAL(DP) :: C1(3),C2(3),PSIG(0:NSPZ+10),
     .          BUFFER(NCHOR,NCHEN),ESTART(NCHOR),ENDFIT(NCHOR),
     .          FP(6), DUM(9)
      REAL(DP) :: ZE1, ZE2, ZSCALE, ZZ, EIRENE_SLOPE, STEIG, PMI, PMA,
     .            XMI, XMAX, XMIN, ZSI, TIMAX, ZE, SUMM, ADD, FAC32,
     .            TEF, DEF, RCMIN, RCMAX, DE, TE, ZDS, RATE, CHKSUM
     .           ,summt,addt, XMA
      REAL(DP) :: EIRENE_FTABRC1
      INTEGER :: I1, I2, IN, I, IS, NAC2, NBC2, ICHRD, IPVOT, NCHNI,
     .           ISK, JSK, IFIRST, JEN, NSPI, ISTR, ICOUNT, KK, IR,
     .           KREC, IRRC, MAXREC, IFLAG, IPLOTS, ILTXT, ISPC,
     .           ICELL, JFEXMN, JFEXMX
      LOGICAL :: NLVL(0:NSTRAI),LCHOR
      CHARACTER(48) :: TX(14)
      CHARACTER(8) :: FILNAM
      CHARACTER(4) :: H123
      CHARACTER(9) :: REAC
      CHARACTER(3) :: CRC
      TYPE(CELL_INFO), POINTER :: FIRST, CUR
C
      ISTRA=IISTR
      NCHNI=IABS(NCHENI)
      IF (NCHTAL(ICHORI).EQ.2) NCHNI=1
C
 
C
      IF (NCHTAL(ICHORI) .NE. 3) THEN
      IF (ISTRA.EQ.IESTR) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=ISTRA
        CALL EIRENE_RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL EIRENE_SYMET(ESTIMV,NTALV,NRAD,NR1ST,NP2ND,NT3RD,
     .               NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
        IESTR=ISTRA
        CALL EIRENE_RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL EIRENE_SYMET(ESTIMV,NTALV,NRAD,NR1ST,NP2ND,NT3RD,
     .               NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
      ELSE
        WRITE (iunout,*) 'ERROR IN DIAGNO: DATA FOR STRATUM ISTRA= ',
     .                   ISTRA
        WRITE (iunout,*)
     .    'ARE NOT AVAILABLE. LINE INTEGRATION ABANDONNED'
        RETURN
      ENDIF
      END IF
 
 
!  EVALUATE SPECTRA COLLECTED IN THE CELLS ALONG THE LINE OF SIGHT
 
      IF ((NCHTAL(ICHORI) == 4) .AND. NLSTCHR(ICHORI)) THEN
 
        if (.not.associated(traj(ichori)%trj%cells)) then
           write (iunout,*) ' HAUE !! '
           return
        end if
        first => traj(ichori)%trj%cells
        cur => first
 
C  RADIATIVE TRANSITION RATES (1/S)
C  BALMER ALPHA
        FAC32=4.410E7
 
        FILNAM='AMJUEL  '
        H123='H.12'
        CRC='OT '
        FP = 0._DP
        RCMIN = -HUGE(1._DP)
        RCMAX =  HUGE(1._DP)
        JFEXMN = 0
        JFEXMX = 0
C
C  H(n=3)/H(n=1)
        REAC='2.1.5a   '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL EIRENE_SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
 
        chksum = 0._dp
        do
          ispc = cur%no_spect
          icell = cur%no_cell
          zds = cur%flight
 
          IF (NSTGRD(ICELL) > 0) GOTO 500
 
          IF (LGVAC(ICELL,NPLS+1)) GOTO 500
 
          TE=TEIN(ICELL)
          DE=DEIN(ICELL)
 
          DEF=LOG(DE*1.D-8)
          TEF=LOG(TE)
 
          call EIRENE_dbl_poly(REACDAT(NREACI+1)%OTH%POLY%DBLPOL,
     .                  TEF, DEF, RATE, DUM, 1, 9,
     .                  RCMIN, RCMAX, FP, JFEXMN, JFEXMX)
          rate = exp(rate)
 
          fuffer(ichori,1:ncheni) = fuffer(ichori,1:ncheni) +
     .                              zds * rate * FAC32 /(4.*PIA) *
     .                              estiml(ispc)%pspc%spc(1:ncheni)
 
          chksum = chksum + zds * rate * FAC32 /(4.*PIA) *
     .                      estiml(ispc)%pspc%spcint
 
 500      continue
          cur => cur%nextc
!pb associated with two arguments tests if both arguments point to the same target
          if (associated(cur,first)) exit
        end do
 
        write (iunout,*) ' chord ', ichori
        write (iunout,*) ' integral along chord calculated from '
        write (iunout,*) ' spectra integrals in the cells '
        write (iunout,*) ' integral ',chksum
        return
 
      END IF
C
C  PREPARE DIRECT (PRIMARY) EMISSIVITY FROM SOURCE, INTO LINE OF SIGHT
C  (SCATTERED (SECONDARY) CONTRIBUTION WILL BE DONE IN SUBR. SIGCX, SIGRAD, SIGH
 
      NLVL=.FALSE.
      DO 10 ISTR=1,NSTRAI
C  PRIMARY SOURCE CONTRIBUTION, REFER TO INPUT BLOCK 7
C  HOWEVER, THIS STRATUM NEED NOT NECESSARILY HAVE TO ACTIVE
        NLVL(ISTR)=NLVOL(ISTR).AND.NLPLS(ISTR)
10    CONTINUE
      NLVL(0)=ANY(NLVL(1:NSTRAI))
C
      IF (NCHTAL(ICHORI).EQ.1) THEN
C  FOR CX SIGNAL:  TO BE WRITTEN
C     CALL ZEROA2(RECADD,NATM,NRAD)
C     IF (NLVL) THEN
C       WRITE (iunout,*) 'WARNING:'
C       WRITE (iunout,*) 'VOLUME RECOMBINATION CONTRIBUTION TO SIGNAL:'
C       WRITE (iunout,*) 'FIRST GENERATION CONTRIBUTION: TO BE WRITTEN'
C  RECADD: #/S/CM**3
C  RECADD*ELCHA*VOL: AMP/CELL
C       DO 100 IPLS=1,NPLSI
C         DO 100 KREC=1,NPRCI(IPLS)
C           IATM=NATPRC(KREC)
C           IF (IATM.LE.0.OR.IATM.GT.NATMI) GOTO 100
C           DO 101 IR=1,NSBOX
C             RECADD(IATM,IR)=RECADD(IATM,IR)+
C    .                        TABRC1(KREC,IR)*DIIN(IPLS,IR)*ELCHA
101         CONTINUE
100     CONTINUE
        write (iunout,*) 'sgnal, cx: ichord,istra,sum ',
     .                      ichori,istra
        write (iunout,*) 'volumetric emission to be written'
C     ENDIF
 
C  FOR RADIANCE OF LINE ISP=IPHOT, IN STRATUM ISTR
      ELSEIF (NCHTAL(ICHORI).EQ.3) THEN
        MAXREC=SUM(NPRCI(1:NPLSI))
        ALLOCATE(RECADD(MAXREC,NRAD))
        ALLOCATE(INTADD(3,MAXREC))
        RECADD=0._DP
        INTADD=0
        ICOUNT=0
        NCTSIG=0
        IF (NLVL(ISTRA)) THEN
C  RECADD: #/S/CM**3
C  RECADD*ELCHA*VOL: AMP/CELL
          IF (ISTRA.EQ.0) THEN
            WRITE (iunout,*) 'ISTRA=0 NOT READY IN SUBR. SGNAL '
            DEALLOCATE(RECADD)
            DEALLOCATE(INTADD)
            LCHOR=.FALSE.
            RETURN
          ENDIF
          IF (NSPEZ(ISTRA).LE.0.OR.NSPEZ(ISTRA).GT.NPLSI) THEN
            WRITE (iunout,*) 'NSPEZ(ISTRA) OUT OF RANGE IN SUBR. SGNAL '
            WRITE (iunout,*) 'ISTRA, NSPEZ ', ISTRA, NSPEZ(ISTRA)
            DEALLOCATE(RECADD)
            DEALLOCATE(INTADD)
            LCHOR=.FALSE.
            RETURN
          ENDIF
          IPLS=NSPEZ(ISTRA)
          IFLAG=0
          DO 130 KREC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,KREC)
            IPHOT=NPHPRC(IRRC)
            IF (IPHOT.EQ.0) IPHOT=NPHPRC_2(IRRC)
            IF (IPHOT.LE.0.OR.IPHOT.GT.NPHOTI) GOTO 130
            IF (IPHOT.NE.ISP) THEN
              IF (NPRCI(IPLS).EQ.1) THEN
                WRITE (IUNOUT,*) 'IN SGNAL: IPHOT.NE.ISP '
                WRITE (IUNOUT,*) 'ICHORD, IRRC, IPHOT, ISP ',
     .                            ICHORI, IRRC,IPHOT,ISP
                WRITE (IUNOUT,*) 'EMISSION FROM THIS LINE IS SKIPPED'
              ELSE
C  THIS PARTICULAR IRRC DOES NOT FIT TO SPECIES ISP FOR PRESENT CHORD
C  SEARCH FOR NEXT IRRC FOR THIS SAME SOURCE PARTICLE IPLS=NSPEZ(ISTRA)
              ENDIF
              GOTO 130
            ENDIF
            KK=NREARC(IRRC)
            ICOUNT=ICOUNT+1
            INTADD(1,ICOUNT)=KK
            INTADD(2,ICOUNT)=IPLS
            INTADD(3,ICOUNT)=IPHOT
            SUMM=0._DP
            summt=0.
            DO 131 IR=1,NSBOX
              ADD=0._DP
              ADDt=0._DP
              IF (NSTGRD(IR).EQ.0.AND..NOT.LGVAC(IR,IPLS)) THEN
                IF (NSTORDR >= NRAD) THEN
                  ADD=TABRC1(IRRC,IR)*DIIN(IPLS,IR)
                ELSE
                  ADD=EIRENE_FTABRC1(IRRC,IR)*DIIN(IPLS,IR)
                END IF
              else
                IF (NSTORDR >= NRAD) THEN
                  ADDt=TABRC1(IRRC,IR)*DIIN(IPLS,IR)
                ELSE
                  ADDt=EIRENE_FTABRC1(IRRC,IR)*DIIN(IPLS,IR)
                END IF
              END IF
              RECADD(ICOUNT,IR)=ADD
              SUMM=SUMM+ADD*VOL(IR)*ELCHA
              SUMMt=SUMMt+ADDt*VOL(IR)*ELCHA
131         CONTINUE
            IFLAG=1
            write (iunout,*) 'sgnal, rad: ichord,istra,icount,',
     .                                                     'vol-source',
     .                            ichori,istra,icount,summ,summt
            write (iunout,*) 'ipls, irrc ',texts(nspami+ipls),' ',irrc
            write (iunout,*) 'iphot, isp ',texts(iphot),' ',texts(isp)
130       CONTINUE
          IF (IFLAG.EQ.0) THEN
C  NO EMISSION FOUND FOR THIS STRATUM, turn off this chord
            LCHOR=.FALSE.
            write (iunout,*) 'sgnal, rad: ichord, istra, is turned off',
     .                               ichori, istra
            write (iunout,*) 'no photon emission rate found '
            DEALLOCATE(RECADD)
            DEALLOCATE(INTADD)
            return
          ENDIF
          NCTSIG=ICOUNT
        ELSE
c  nlvl is not true:
          LCHOR=.FALSE.
          write (iunout,*) 'sgnal, rad: ichord, istra, is turned off',
     .                               ichori, istra
          write (iunout,*) 'specified stratum is not a volume source '
          DEALLOCATE(RECADD)
          DEALLOCATE(INTADD)
          return
        ENDIF
      ELSEIF (NCHTAL(ICHORI).EQ.2) THEN
        write (iunout,*) 'sgnal, emis: ichord,istra ',
     .                            ichori,istra
        write (iunout,*) 'volumetric line emission '
        write (iunout,*) 'contribution no. isp ',isp
      ENDIF
C
C     CALCULATE SIGNAL STRENGTHS
C
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
C
C  ENERGY LOOP (IF ANY)
C
      PMA=-1.E30
      PMI=1.E30
      XMA=-1.E30
      XMI=1.E30
      IF (NCHTAL(ICHORI).EQ.1) NSPI=NATMI
      IF (NCHTAL(ICHORI).EQ.2) NSPI=10
      IF (NCHTAL(ICHORI).EQ.3) NSPI=NPHOTI
      IF (NCHTAL(ICHORI).EQ.10) NSPI=NSPZ
      PSIG = 0._DP
      IFIRST=0
      DO 231 JEN=1,NCHNI
        ZE=ENERGY(JEN)
        CALL EIRENE_LININT
     .  (IFIRST,ICHORI,C1,C2,ICHRD,IPVOT,NBC2,NAC2,ZE,
     .               PSIG,TIMAX,ISP,NSPI,JEN,NCHNI)
        IFIRST=1
        IF (ISP.GT.0.AND.ISP.LE.NSPI) THEN
C  SINGLE SPECIES INDEX ISP
          BUFFER(ICHORI,JEN)=PSIG(ISP)
        ELSEIF (ISP.EQ.0) THEN
C  SUM OVER SPECIES INDEX
          ZSI=0.
          DO 239 IS=1,NSPI
            ZSI=ZSI+PSIG(IS)
239       CONTINUE
          BUFFER(ICHORI,JEN)=ZSI
          WRITE (80,'(I6,3ES12.4)') ICHORI,C2
          WRITE (80,'(6ES12.4)') PSIG(0:NSPI)
        ELSE
          WRITE (iunout,*) 'ERROR IN SUBR. SGNAL: ISP= ',ISP
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
C  PROCESS DATA FROM LINE INTEGRAL ROUTINES INTO REQUESTED DATA & UNITS
C
        IF (NCHTAL(ICHORI).EQ.1) THEN
C  LINE INTEGRAL: CX ATOMS/SEC/CM**2/EV/STERAD
C  THE NUMERICAL FACTOR 1./11.137 ARISES FROM A TRANSFORMATION
C  OF A MAXWELLIAN VELOCITY DISTRIBUTION TO A MAXW. ENERGY DISTR.
C  1./11.137=0.5*(1./PI)**1.5, IN SIGCX
          FUFFER(ICHORI,JEN)=BUFFER(ICHORI,JEN)/11.137
        ELSEIF (NCHTAL(ICHORI).EQ.2) THEN
C  LINE INTEGRAL: PHOTONS/SEC/CM**2/STERAD (EMISSIVITY)
          FUFFER(ICHORI,JEN)=BUFFER(ICHORI,JEN)/(4.*PIA)
        ELSEIF (NCHTAL(ICHORI).EQ.3) THEN
C  LINE INTEGRAL: PHOTONS/SEC/CM**2/EV/STERAD (SPECTRAL RADIANCE)
          FUFFER(ICHORI,JEN)=BUFFER(ICHORI,JEN)/(4.*PIA)
        ELSEIF (NCHTAL(ICHORI).EQ.10) THEN
C  LINE INTEGRAL: USER SUPPLIED INTEGRAND ALONG LINE OF SIGHT
          FUFFER(ICHORI,JEN)=BUFFER(ICHORI,JEN)
        ENDIF
231   CONTINUE
C
C  ENERGY LOOP FINISHED
C
C
C  FIT A STRAIGHT LINE TO SPECTRUM, BETWEEN ESTART AND ENDFIT
C  TO BE DONE ONLY FOR NCHTAL=1
C
      IF ((NCHTAL(ICHORI).NE.1) .OR. (NCHENI == 1)) GOTO 300
 
      ESTART(ICHORI)=NSPINI(ICHORI)*TIMAX
      ENDFIT(ICHORI)=NSPEND(ICHORI)*TIMAX
      TINP(ICHORI)=TIMAX
      IF (TRCSIG) THEN
        WRITE (iunout,*) 'FITTING RANGE: E1--E2, TIMAX'
        WRITE (iunout,*) ESTART(ICHORI),'--',ENDFIT(ICHORI),'  ',TIMAX
      ENDIF
C
C  FIND MAX. VALUE OF SIGNAL
      CALL EIRENE_MAXMN2(BUFFER,NCHOR,ICHORI,ICHORI,1,NCHNI,XMIN,XMAX)
C  SCALE RESULT
      IF (XMAX.GT.0.) GOTO 235
      CALL EIRENE_MASAGE
     .  ('NO SLOPE IN SIGNAL, BECAUSE MAX(BUFFER).LE.0   ')
      PLSPEC=.FALSE.
      RETURN
235   ZSCALE=1./XMAX
      DO 233 I=1,NCHNI
        BUFFER(ICHORI,I)=BUFFER(ICHORI,I)*ZSCALE
        ZZ=MAX(1.E-10_DP,BUFFER(ICHORI,I))
233     BUFFER(ICHORI,I)=LOG(ZZ)
C
C  CURVE FITTING
C
C  MIN. ENERGY FOR FIT
      ZE1=ESTART(ICHORI)
C  MAX. ENERGY
      ZE2=ENDFIT(ICHORI)
C
C  FIND ELEMENTS
      DO 241 JEN=2,NCHNI
        I1=JEN
        IF (ENERGY(I1).GE.ZE1) GO TO 242
241   CONTINUE
242   CONTINUE
C
      DO 243 JEN=I1,NCHNI
        I2=JEN
        IF (ENERGY(I2).GE.ZE2) GO TO 244
243   CONTINUE
244   CONTINUE
C
C   NUMBER OF POINTS FOR FITTING
      IN=I2-I1+1
      IF (IN.LT.2) THEN
        CALL EIRENE_MASAGE
     .  ('WRONG ENERGY RANGE FOR CURVE FITTING IN SIGNAL ')
        CALL EIRENE_MASJ1 ('CHORD.NO. I=           ',ICHORI)
        WRITE (iunout,*) 'I1,I2,IN ',I1,I2,IN
        TILINE(ICHORI)=0.
      ELSE
C  FITTED RESULT
        STEIG=EIRENE_SLOPE(IN,ICHORI,I1,ENERGY,BUFFER,NCHOR,NCHNI)
        IF (STEIG.GE.0.) THEN
          TILINE(ICHORI)=0.
        ELSE
          TILINE(ICHORI)=-1./STEIG
        ENDIF
      ENDIF
C
      IF (TRCSIG) THEN
        WRITE (iunout,*) 'ICHORI,TILINE(ICHORI) ',ICHORI,TILINE(ICHORI)
      ENDIF
 
300   CONTINUE
C
      IF (NCHTAL(ICHORI).EQ.3) THEN
        DEALLOCATE(RECADD)
        DEALLOCATE(INTADD)
      ENDIF
      RETURN
      END
