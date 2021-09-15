C
C
      SUBROUTINE SIGCX(IFIRST,JJJ,ZDS,PEN,PSIG,TIMAX,ARGST)
C
C  INPUT:
C          IFIRST: FLAG FOR INITIALISATION
C          NCELL:  INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C          PEN:    ENERGY (EV) AT WHICH CX FLUXES ARE TO BE EVALUATED
C  OUTPUT: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO:
C          THE CX FLUX PSIG(IATM),IATM=0,NATMI OF ENERGY PEN (EV)
C          THE MAX. ION TEMP. TIMAX ALONG LINE OF SIGHT
C          THE INTEGRANT ARGST SUCH THAT INTEGR.(ARGST*DL) = PSIG
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE COMSIG
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CGEOM
      USE COMPRT
      USE CLGIN
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IFIRST, JJJ
      REAL(DP), INTENT(IN) :: ZDS, PEN
      REAL(DP), INTENT(IN OUT) :: PSIG(0:NSPZ+10), TIMAX
      REAL(DP), INTENT(IN OUT) :: ARGST(0:NSPZ+10,NRAD)
      REAL(DP), ALLOCATABLE, SAVE :: ZARG2(:), ZARG3(:)
      REAL(DP) :: ZLAMB(0:NATM),
     .          SIGTTT(0:NATM,NPLS), CFLAG(7,3), ZNI(NRCX)
      REAL(DP) :: ELAB, CXS, CXRATE, CROSS, PVELQ0, HEB, VREL,
     .          VRELQ, ZEXP2, ATTENU, ZEXP3, SIGADD, ARGU, ZEXP1,
     .          ZMAX, VXS, VYS, VZS, ZTII, ZARG1, ZTI32, FPATHA, ZTI,
     .          TTARG, EDR
      INTEGER :: IREAC, IAT, NCELC, ICELL, IRCX, IACX, IPLSTI
      INTEGER :: ISPC, IE
C
      DATA ZMAX/40./

      IF (IFIRST.EQ.0) THEN
        ALLOCATE (ZARG2(0:NATM))
        ALLOCATE (ZARG3(0:NATM))
        TIMAX=0.
        DO 100 IATM=0,NATMI
          PSIG(IATM)=0.
          ZARG2(IATM)=0.
          ZARG3(IATM)=0.
          DO 101 ICELL=1,NSBOX
            ARGST(IATM,ICELL)=0.
101       CONTINUE
100     CONTINUE
      ELSEIF (IFIRST.EQ.1) THEN
C
C  MEAN FREE PATH LENGTH FOR ATTENUATION FACTOR
C  IONS IPLS WITH  (SHIFTED, VXIN,VYIN,VZIN) MAXWELLIAN
C  AT KT=TIIN(IPLS,NCELL),
C  ALL SPECIES IN PLASMA
C  ELECTRONS WITH  MAXWELLIAN AT KT=TEIN(NCELL),
C  MONOENERGETIC NEUTRAL BEAM, VELOCITY VEL (CM/SEC), SPECIES IATM
C  WITH SPEED UNIT VECTOR (-VELX,-VELY,-VELZ)
C
        NCELC=NCLTAL(NCELL)

c       DO 200 IATM=1,NATMI
c         VEL=SQRT(PEN/RMASSA(IATM))*CVELAA
c         VXS=VELX
c         VYS=VELY
c         VZS=VELZ
c         VELX=-VELX
c         VELY=-VELY
c         VELZ=-VELZ
c         ZLAMB(IATM)=FPATHA(NCELL,CFLAG,1,1)
c         VELX=VXS
c         VELY=VYS
c         VELZ=VZS
200     CONTINUE
C
C  CX-COLLISION FREQUENCY (1/SEC)
C  MONOENERGETIC ION BEAM, ENERGY PEN (EV), SPEED UNIT VECTOR:
C  -(VELX,VELY,VELZ) AND SPECIES IPLS=1,NPLSI
C  NEUTRALS WITH MAXWELLIAN KT=2/3*EDENA/PDENA, SPECIES IATM=1,NATMI
C  ASSUME: POST COLLISION NEUTRAL = PRE COLLISION ION: PEN, -(VELX,VELY,VELZ)
C
c       DO 300 IATM=0,NATMI
c         DO 300 IPLS=1,NPLSI
c           SIGTTT(IATM,IPLS)=0.
300     CONTINUE
        ZNI = 0._DP
C
        DO 311 IATM=1,NATMI
          IF (LGACX(IATM,0,0).EQ.0) GOTO 311
          IF (LGVAC(NCELL,0).OR.PDENA(IATM,NCELC).LE.0.) GOTO 311

c  find value of spectrum of iatm, in cell ncell, with energy pen
          if (lspccll(ncell)) then
            zni(1) = 0._dp
!            do ispc=1,nadspc
!              if ((estiml(ispc)%pspc%isrfcll==2) .and.
!     .            (estiml(ispc)%pspc%ispcsrf==ncell) .and.
!     .            (estiml(ispc)%pspc%iprtyp==1) .and.
!     .            (estiml(ispc)%pspc%iprsp==iatm)) then
!              IF (PEN < ESTIML(ISPC)%PSPC%SPCMIN) THEN
!                IE = 0
!              ELSEIF (PEN >= ESTIML(ISPC)%PSPC%SPCMAX) THEN
!                IE = ESTIML(ISPC)%PSPC%NSPC + 1
!              ELSE
!                IE = (PEN - ESTIML(ISPC)%PSPC%SPCMIN) *
!     .                ESTIML(ISPC)%PSPC%SPCDELI + 1
!              END IF
!              zni(1)=estiml(ispc)%pspc%spc(ie)
!              exit
!              end if
!            enddo
c  hard wired: use background spectrum for ipls=5 (d_3g)
            do ispc=1,nback_spec
              if ((back_spec(ispc)%pspc%isrfcll==2) .and.
     .            (back_spec(ispc)%pspc%ispcsrf==ncell) .and.
     .            (back_spec(ispc)%pspc%iprtyp==4) .and.
     .            (back_spec(ispc)%pspc%iprsp==5)) then
              IF (PEN < back_spec(ISPC)%PSPC%SPCMIN) THEN
                IE = 0
              ELSEIF (PEN >= back_spec(ISPC)%PSPC%SPCMAX) THEN
                IE = back_spec(ISPC)%PSPC%NSPC + 1
              ELSE
                IE = (PEN - back_spec(ISPC)%PSPC%SPCMIN) *
     .                back_spec(ISPC)%PSPC%SPCDELI + 1
              END IF
              zni(1)=back_spec(ispc)%pspc%spc(ie)
              exit
              end if
            enddo
          else
            write (6,*)  'kein vorabgedoens '
            call exit_own(1)
          endif
C  NEUTRAL DRIFT VELOCITIES (XDR,YDR,ZDR) NOT YET ACCOUNTED FOR.
C  INCIDENT ION BEAM (PEN) MUST BE CONVERTED INTO A FRAME IN WHICH
C  NEUTRAL DRIFT=0. NEXT LINE MUST CONTAIN NEUTRAL DRIFT ENERGY.
c         EDR=0.
c         TTARG=2./3.*(EDENA(IATM,NCELC)-EDR)/PDENA(IATM,NCELC)
c         ZTI=TTARG*8./PIA*CVEL2A*CVEL2A
c         DO 310 IACX=1,NACXI(IATM)
c           IRCX=LGACX(IATM,IACX,0)
c           IPLS=LGACX(IATM,IACX,1)
C           IAT_IN = IATM
C           IPL_IN = IPLS
C           IAT_OUT= N1STX(IRCX,2)
C           IPL_OUT= N2NDX(IRCX,2)
C
C  LOCAL ION TEMPERATURE
C
c           IPLSTI = MPLSTI(IPLS)
c           ZTII=TIIN(IPLSTI,NCELL)
c           TIMAX=MAX(TIMAX,ZTII)
C
C  ZNI=PROBABLILITY FOR EMISSION OF ATOMS IAT_OUT WITH E=PEN:
C          FROM MAXW. ENERGY DISTR. OF IONS IPL_IN GOING INTO CX
C  ZNI: 1/EV
c           ZTI32=1.0/(ZTII*SQRT(ZTII))
c           ZARG1=PEN/ZTII
c           ZEXP1=0.
c           IF(ZARG1.LE.ZMAX) ZEXP1=EXP(-ZARG1)
c           ZNI(IRCX)=(ZEXP1*ZTI32)*SQRT(PEN)
C  cx process data for ircx assume: beam iatm, ions: ipls
c           VEL=SQRT(PEN/RMASSA(IATM))*CVELAA
c           PVELQ0=VEL*VEL
c           HEB=LOG10(PVELQ0*CVELI2)
C
c           VRELQ=ZTI/RMASSP(IPLS)+PVELQ0
c           VREL=SQRT(VRELQ)
C  STRICTLY: A BEAM (ION, IPL_IN) MAXWELLIAN (NEUTRAL, IAT_OUT)
C  RATE COEFFICIENT SHOULD BE COMPUTED
C  this should be obtained from <sigma v>(IRCX) with proper
c  choice of Ti_scaled and E_scaled
C  USE APPROXIMATION: <SIGMA V> APPROX SIGMA(V_EFF)*V_EFF
C  IS CX CROSS SECTION AVAILABLE ?
c           IREAC=MODCOL(3,1,IRCX)
c           IF (IREAC.EQ.-1.OR.IREAC.GT.0) THEN
c             ELAB=LOG(VRELQ)+DEFCX(IRCX)
c
c             CXS=CROSS(ELAB,IREAC,IRCX,'SIGCX ')
c             CXRATE=CXS*VREL
C
c             SIGTTT(IATM,IPLS)=CXRATE*PDENA(IATM,NCELC)
c           ELSE
c             WRITE (iunout,*) 'CROSS SECTION NOT AVAILABLE IN SIGCX  '
c           ENDIF
310       CONTINUE
C
311     CONTINUE
C
C  UP TO THIS POINT, IATM IS THE SPECIES INDEX OF THE ATOM GOING INTO
C  CX WITH BULK ION IPLS, AND SIGTTT(IATM,IPLS) IS THE CORRESPONDING RATE.
C  FOR THE SPECTRUM, NOW, THE ATOM SPECIES AFTER CX MATTERS 
C  (E.G., IF : H + D+ --> D + H+
C  THEN, UP TO NOW, IATM WAS H.  FROM NOW ON IT MUST STAND FOR D.)
C
C  ATTENUATION FACTOR:
C  FOR ALL CX ATOMS WITH SPECIES INDEX IATM, TRAVELLING IN IPLS
        DO 400 IATM=1,NATMI
c         ZARG3(IATM)=ZARG3(IATM)+ZARG2(IATM)
c         ZARG2(IATM)=ZDS/ZLAMB(IATM)
c         ZEXP2=0.0
c         ZEXP3=0.0
c         IF(ZARG2(IATM).LT.ZMAX) ZEXP2=EXP(-ZARG2(IATM))
c         IF(ZARG3(IATM).LT.ZMAX) ZEXP3=EXP(-ZARG3(IATM))
c         ATTENU=ZLAMB(IATM)*ZEXP3*(1-ZEXP2)
C
C  SOURCE-TERM FOR ATOMS IATM
C
C  1.) CX OF IONS IPLS WITH ATOMS IAT_IN, RESULTING IN IATM=IAT_OUT
C  2.) DIRECT  FROM PRIMARY SOURCE (E.G. RECADD: RECOMBINATION)
C  3.) DIRECT  FROM NON-CX  SECUNDARY SOURCE (E.G.  WALL REFLECTION INTO PEN, LINE)
C
c       SIGADD=0.
C  CONTRIBUTION 1:
c       DO IRCX = 1, NRCX
c         IF ( N1STX(IRCX,2) /= IATM ) CYCLE
c           IF (ANY(LGACX(:,:,0)==IRCX)) THEN
C  IAT IS THAT ATOM SPECIES, WHICH RESULTS IN IATM AFTER CX WITH IPLS
C  FIND PRE-COLLISION ATOM INDEX IAT:
c             DO IAT=1,NATMI
c               DO IACX=1,NACXI(IAT)
c                 IF (LGACX(IAT,IACX,0) == IRCX) THEN
c                   IPLS = LGACX(IAT,IACX,1)
c                   SIGADD=SIGTTT(IAT,IPLS)*DIIN(IPLS,NCELL)+SIGADD
c                 END IF
c             END DO
c             END DO
c           END IF
cdr         ARGU=ZNI(IRCX)*SIGADD
            ARGU=ZNI(1)
            PSIG(IATM)=PSIG(IATM)+ARGU*ZDS
            ARGST(IATM,JJJ)=ARGU
c       END DO
C  CONTRIBUTION 2: TO BE WRITTEN
C         SIGADD=RECADD(IATM,NCELL)  ......
C  CONTRIBUTION 3: TO BE WRITTEN
C         SIGADD= .............
400     CONTINUE
C
      ELSEIF (IFIRST.EQ.2) THEN
        DEALLOCATE (ZARG2)
        DEALLOCATE (ZARG3)
      ENDIF
      RETURN
      END
