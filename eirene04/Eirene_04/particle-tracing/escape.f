CVK 25.02.04: splitting for sputtering is re-introduced, (from eirene_02)
CVK 25.02.04: the ILIIN=-3 option support is re-introduced
CVK 25.02.04: spttot updated (total sputtered flux tally)
CDR 25.02.04: return, return1 for reflected flux tallies moved after call
CDR 25.02.04:                 to upsusr, calc_spectrum (from eirene_02)
C
      SUBROUTINE ESCAPE(PR,SG,*,*,*)
C
C  PROCESS ESCAPING PARTICLES
C  INPUT:
C        LGPART=.TRUE.  UPDATE TALLIES FOR INCIDENT PARTICLES, THEN SURFACE MODEL
C        LGPART=.FALSE. UPDATE TALLIES FOR INCIDENT PARTICLES, THEN STOP.
C
C  RETURN  : STOP TRACK OF THIS PARTICLE TYPE. RETURN TO SUBR. MCARLO
C  RETURN 1: START NEW TRACK OF SAME TYPE IN CALLING PROGRAM FOLNEUT
C            OR FOLION
C  RETURN 2: CONTINUE THIS TRACK IN CALLING PROGRAM (TRANSP. SURFACE)
C  RETURN 3: RESTORE COLLISION DATA, CONDITIONAL EXPECTATION WAS USED
C
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CLOGAU
      USE CRAND
      USE CINIT
      USE CGRID
      USE CSPEZ
      USE CZT1
      USE COMPRT
      USE COMSPL
      USE CLGIN
      USE COUTAU
      USE CTRIG

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PR, SG
      REAL(DP) :: DIWL(NPLS), VPWL(NPLS)
      REAL(DP) :: VXSPTC, VYSPTC, VZSPTC, VSPTC, ESPTC, VXSPTP,
     .          VYSPTP, VZSPTP, ESPTP, VELXS, VELYS, VELZS, VSPTP, TW,
     .          E0TERM, ZEP1, FR2, COSI2, ZVZ, WABS,
     .          CUR, GAMMA, TEWL, VX, VY, VZ, FCHAR, WPR, FMASS,
     .          FLX, YIELD1, YIELD2, VELS, WEIGHS, E0S, ESHET,
     .          VSHETQ, V, VELSH, VC, SHEATH, SPLFLG
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: ISG, ISPZS, I, J, IDIM, MS, IC, IP, ITOLD, ISTS,
     .           ISSPTP, ISSPTC, IPV
      LOGICAL :: NLSPUT, LTRANS
C
C
C  .................
C  .               .
C  .  PERIODICITY: .
C  .................
C
C  CURRENTLY: ONLY IN CASE LEVGEO=1, NLTRZ. HENCE: VEL_OLD=VEL_NEW
C  TO BE WRITTEN: TOROIDICITY AS SPECIAL CASE OF PERIODICITY
      IF (ILIIN(MSURF).GE.4) THEN
C  CONDITIONAL EXPECTATION ESTIMATOR: HAS THIS PARTICLE COLLIDED IN THE VOLUME,
C  BEFORE IT HIT THE WALL?
        IF (ICOL.EQ.1) RETURN 3
C  NO
        IF (.NOT.LGPART) THEN
          WRITE (6,*) 'ERROR AT PERIODICITY SURFACE, LGPART=FALSE '
          RETURN
        ENDIF
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,11)
        RETURN 1
      ENDIF
C
C  .............................
C  .                           .
C  .  UPDATE INCIDENT FLUXES:  .
C  .............................
C
      ITOLD=ITYP
      WPR=WEIGHT*PR
C
C  UPDATE PARTICLE EFFLUX  ONTO SURFACE MSURF
C  UPDATE ENERGY FLUX ONTO SURFACE MSURF
C
C  SPATIAL RESOLUTION ON NON DEFAULT STANDARD SURFACE?
      IF (MSURF.GT.NLIM.AND.NLMPGS.GT.NLIMPS) THEN
        IF (LEVGEO.LE.3) THEN
          ISTS=MSURF-NLIM
          MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
          FLX=FLXOUT(MSURFG)
        ELSE IF (LEVGEO.EQ.4) THEN
          MSURFG=NLIM+NSTS+INSPAT(IPOLGN,MRSURF)
          FLX=FLXOUT(MSURFG)
        ELSE
          MSURFG=0
          FLX=FLXOUT(MSURF)
        END IF
      ELSEIF (MSURF.GT.0) THEN
        MSURFG=0
        FLX=FLXOUT(MSURF)
      ELSE
        WRITE (6,*) 'ERROR IN ESCAPE: MSURF=0. KILL PARTICLE '
        RETURN
      ENDIF
C
C  FOR NON-TRANSPARENT SURFACES: INCIDENT FLUX
C  FOR     TRANSPARENT SURFACES: ONE SIDED FLUX, POSITIVE COMPONENT
C
      IF ((ILIIN(MSURF).LT.0).AND.(SG.LT.0.D0).AND.(ILIIN(MSURF).NE.-3)) 
     .GOTO 10
C
      IF (ITYP.EQ.1) THEN
        IF (LEOTAT) EOTAT(IATM,MSURF)=EOTAT(IATM,MSURF)+E0*WPR
        IF (LPOTAT) POTAT(IATM,MSURF)=POTAT(IATM,MSURF)+WPR
        FMASS=DBLE(NMASSA(IATM))
        FCHAR=DBLE(NCHARA(IATM))
      ELSEIF (ITYP.EQ.2) THEN
        IF (LEOTML) EOTML(IMOL,MSURF)=EOTML(IMOL,MSURF)+E0*WPR
        IF (LPOTML) POTML(IMOL,MSURF)=POTML(IMOL,MSURF)+WPR
        FMASS=DBLE(NMASSM(IMOL))
        FCHAR=DBLE(NCHARM(IMOL))
      ELSEIF (ITYP.EQ.3) THEN
        IF (ILIIN(MSURF).GT.0) THEN
          ESHET=0.D0
C  ACCOUNT FOR ELECTROSTATIC SHEATH AT SURFACE FOR TEST IONS
          IF (FSHEAT(MSURF).LE.0.D0) THEN
            GAMMA=0.
            CUR=0.
            IC=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
            IF (.NOT.LGVAC(IC,NPLS+1)) THEN
              TEWL=TEIN(IC)
              DO 30 IP=1,NPLSI
                IPV=MPLSV(IP)
                IF (INDPRO(4) == 8) THEN
                  CALL VECUSR (2,VX,VY,VZ,IP)
                ELSE
                  VX=VXIN(IPV,IC)
                  VY=VYIN(IPV,IC)
                  VZ=VZIN(IPV,IC)
                ENDIF
                VPWL(IP)=SQRT(VX**2+VY**2+VZ**2)
                DIWL(IP)=DIIN(IP,IC)
30            CONTINUE
              ESHET=NCHRGI(IION)*SHEATH(TEWL,DIWL,VPWL,
     .                                  NCHRGP,GAMMA,CUR,NPLSI,MSURF)
            ENDIF
          ELSE
            IC=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
            IF (.NOT.LGVAC(IC,NPLS+1)) THEN
              TEWL=TEIN(IC)
              ESHET=FSHEAT(MSURF)*TEWL*NCHRGI(IION)
            ENDIF
          ENDIF
C  ADD VELOCITY DUE TO SHEATH ACCELERATION
          IF (ESHET.GT.0.D0) THEN
            VSHETQ=ESHET*RSQDVI(IION)*RSQDVI(IION)
            VC=2.*VEL*(VELX*CRTX+VELY*CRTY+VELZ*CRTZ)
            VELSH=-VC/2.+SQRT(VC*VC/4.+VSHETQ)
            VX=VEL*VELX+VELSH*CRTX
            VY=VEL*VELY+VELSH*CRTY
            VZ=VEL*VELZ+VELSH*CRTZ
            V=SQRT(VX*VX+VY*VY+VZ*VZ)
            VELX=VX/V
            VELY=VY/V
            VELZ=VZ/V
            E0=E0+ESHET
            E0_MEAN=E0
            VEL=SQRT(E0)*RSQDVI(IION)
            EELFI(IION,ISTRA)=EELFI(IION,ISTRA)+ESHET*WPR
          ENDIF
        ENDIF
        IF (LEOTIO) EOTIO(IION,MSURF)=EOTIO(IION,MSURF)+E0*WPR
        IF (LPOTIO) POTIO(IION,MSURF)=POTIO(IION,MSURF)+WPR
        FMASS=DBLE(NMASSI(IION))
        FCHAR=DBLE(NCHARI(IION))
      ELSEIF (ITYP.EQ.0) THEN
        IF (LEOTPHT) EOTPHT(IPHOT,MSURF)=EOTPHT(IPHOT,MSURF)+E0*WPR
        IF (LPOTPHT) POTPHT(IPHOT,MSURF)=POTPHT(IPHOT,MSURF)+WPR
        FMASS=0._dp
        FCHAR=0._dp
      ENDIF
      IF (MSURFG.GT.0) THEN
        IF (ITYP.EQ.1) THEN
          IF (LEOTAT) EOTAT(IATM,MSURFG)=EOTAT(IATM,MSURFG)+E0*WPR
          IF (LPOTAT) POTAT(IATM,MSURFG)=POTAT(IATM,MSURFG)+WPR
        ELSEIF (ITYP.EQ.2) THEN
          IF (LEOTML) EOTML(IMOL,MSURFG)=EOTML(IMOL,MSURFG)+E0*WPR
          IF (LPOTML) POTML(IMOL,MSURFG)=POTML(IMOL,MSURFG)+WPR
        ELSEIF (ITYP.EQ.3) THEN
          IF (LEOTIO) EOTIO(IION,MSURFG)=EOTIO(IION,MSURFG)+E0*WPR
          IF (LPOTIO) POTIO(IION,MSURFG)=POTIO(IION,MSURFG)+WPR
        ELSEIF (ITYP.EQ.0) THEN
          IF (LEOTPHT) EOTPHT(IPHOT,MSURFG)=EOTPHT(IPHOT,MSURFG)+E0*WPR
          IF (LPOTPHT) POTPHT(IPHOT,MSURFG)=POTPHT(IPHOT,MSURFG)+WPR
        ENDIF
      ENDIF
      ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
10    CONTINUE
C
C  ADDITIONAL INCIDENT SURFACE TALLIES
      IF (NADSI.GE.1) CALL UPSUSR (WPR,1)
      IF (NADSPC.GE.1) CALL CALC_SPECTRUM (WPR,1)
C
C  STOP TRAJECTORY, FOR SOME REASON IN SUBROUTINE ADDCOL OR STDCOL
C
      IF (.NOT.LGPART) RETURN
C
C   ........................
C   .                      .
C   .  CALL SPUTTER MODEL  .
C   ........................
C
C
C  TENTATIVELY ASSUME: NO SPUTTERED PARTICLES WILL BE FOLLOWED
      NLSPUT=.FALSE.
C
      WGHTSP=0.
      WGHTSC=0.
      YIELD1=0.
      YIELD2=0.
      ISSPTP=0
      ISSPTC=0
C
      IF (ILSPT(MSURF).NE.0) THEN
C
C  AT THIS POINT: ILIIN.EQ.1
C  SAVE PARAMETERS OF INCIDENT PARTICLE
        E0S=E0
        WEIGHS=WEIGHT
        VELS=VEL
        VELXS=VELX
        VELYS=VELY
        VELZS=VELZ
        ISPZS=ISPZ
        CALL SPUTR1(WMINS,FMASS,FCHAR,FLX,
     .              ISRS(ISPZ,MSURF),
     .              YIELD1,
     .              ISSPTP,ESPTP,VSPTP,VXSPTP,VYSPTP,VZSPTP,
     .              ISRC(ISPZ,MSURF),
     .              YIELD2,
     .              ISSPTC,ESPTC,VSPTC,VXSPTC,VYSPTC,VZSPTC)
        NLSPUT=YIELD1.GT.0..OR.YIELD2.GT.0.
        WGHTSP=WPR*YIELD1
        WGHTSC=WPR*YIELD2
C
C  UPDATE SPUTTERED FLUX IF AVAILABLE. SORTED BY INCIDENT PARTICLE TYPE
C
        IF (NLSPUT) THEN
C  UPDATE TOTAL SPUTTERED FLUX TALLY
          IF (LSPTTOT) SPTTOT(MSURF)=SPTTOT(MSURF)+WGHTSP+WGHTSC
          IF (ITYP.EQ.1) THEN
            IF (LSPTAT) 
     .        SPTAT(IATM,MSURF)=SPTAT(IATM,MSURF)+WGHTSP+WGHTSC
          ELSEIF (ITYP.EQ.2) THEN
            IF (LSPTML) 
     .        SPTML(IMOL,MSURF)=SPTML(IMOL,MSURF)+WGHTSP+WGHTSC
          ELSEIF (ITYP.EQ.3) THEN
            IF (LSPTIO) 
     .       SPTIO(IION,MSURF)=SPTIO(IION,MSURF)+WGHTSP+WGHTSC
C         ELSEIF (ITYP.EQ.4) THEN
C           IF (LSPTPL) 
C     .       SPTPL(IPLS,MSURF)=SPTPL(IPLS,MSURF)+WGHTSP+WGHTSC
          ELSEIF (ITYP.EQ.0) THEN
            IF (LSPTPHT) 
     .        SPTPHT(IPHOT,MSURF)=SPTPHT(IPHOT,MSURF)+WGHTSP+WGHTSC
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LSPTTOT) SPTTOT(MSURFG)=SPTTOT(MSURFG)+WGHTSP+WGHTSC
            IF (ITYP.EQ.1) THEN
              IF (LSPTAT) 
     .          SPTAT(IATM,MSURFG)=SPTAT(IATM,MSURFG)+WGHTSP+WGHTSC
            ELSEIF (ITYP.EQ.2) THEN
              IF (LSPTML) 
     .          SPTML(IMOL,MSURFG)=SPTML(IMOL,MSURFG)+WGHTSP+WGHTSC
            ELSEIF (ITYP.EQ.3) THEN
              IF (LSPTIO) 
     .          SPTIO(IION,MSURFG)=SPTIO(IION,MSURFG)+WGHTSP+WGHTSC
C           ELSEIF (ITYP.EQ.4) THEN
C              IF (LSPTPL) 
C     .         SPTPL(IPLS,MSURFG)=SPTPL(IPLS,MSURFG)+WGHTSP+WGHTSC
            ELSEIF (ITYP.EQ.0) THEN
              IF (LSPTPHT) 
     .          SPTPHT(IPHOT,MSURFG)=SPTPHT(IPHOT,MSURFG)+WGHTSP+WGHTSC
            ENDIF
          ENDIF
        ENDIF
C
      ENDIF
C
C
C   ...................................................................
C   .                                                                 .
C   .  SEMI-TRANSPARENT SURFACE, FOR A FRACTION "TRANSP" OF THE FLUX  .
C   ...................................................................
C
C
      LTRANS=.FALSE.
      IF (TRANSP(ISPZ,1,MSURF).GT.0.D0.OR.
     .    TRANSP(ISPZ,2,MSURF).GT.0.D0) THEN
C
C  AT THIS POINT: ILIIN(MSURF).GT.0
C
        IF (SG.GT.0) ISG=1
        IF (SG.LT.0) ISG=2
        LTRANS=RANF_EIRENE( ).LE.TRANSP(ISPZ,ISG,MSURF)
        IF (LTRANS) THEN
C  A NON TRANSPARENT SURFACE IS MADE TRANSPARENT FOR THIS
C  PARTICULAR PARTICLE
C  STANDARD OR ADDITIONAL SURFACE?
          MS=MSURF
          IF (MSURF.GT.NLIM) THEN
            ISTS=MSURF-NLIM
            MS=-ISTS
            IF (INUMP(ISTS,1).NE.0) IDIM=1
            IF (INUMP(ISTS,2).NE.0) IDIM=2
            IF (INUMP(ISTS,3).NE.0) IDIM=3
C  CELL NUMBER SWITCHES LIKE A TRANSPARENT DEFAULT STANDARD SURFACE IN STDCOL
            IF (IDIM.EQ.1) NRCELL=NRCELL+NINCX
            IF (IDIM.EQ.2) THEN
              NPCELL=NPCELL+NINCY
              IPOLG=MPSURF
            END IF
            IF (IDIM.EQ.3) NTCELL=NTCELL+NINCZ
          ENDIF
          IF (NLTRC) THEN
            CALL LEER(1)
            WRITE (6,*) 'SURFACE MSURF= ',MS,' IS MADE TRANSPARENT'
            WRITE (6,*) 'ORIENTATION, SPECIES: ',ISG,ISPZ
          ENDIF
        ENDIF
C
      ENDIF
C
C
C   ...................................
C   .                                 .
C   .  PERFECTLY ABSORBING SURFACES   .
C   ...................................
C
C
C  NOTHING ELSE TO BE DONE, RETURN
C
      IF (ILIIN(MSURF).EQ.2.AND..NOT.LTRANS) THEN
        IF (LSPUMP) SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WPR
        WEIGHT=0.D0
        LGPART=.FALSE.
        IF (ICOL.EQ.1) RETURN 3
        RETURN
      ENDIF
C
C   ..............................................
C   .                                            .
C   .   MIRROR, OR SEMI-TRANSPARENT SURFACE      .
C   .   REEMITTED FLUX=INCOMING FLUX AND RETURN  .
C   ..............................................
C
      IF (LTRANS.OR.ILIIN(MSURF).EQ.3) THEN
C
C ITOLD=ITNEW=ITYP
C
        IF (ITYP.EQ.1) THEN
          IF (LERFAAT) ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WPR
          IF (LPRFAAT) PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WPR
        ELSEIF (ITYP.EQ.2) THEN
          IF (LERFMML) ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WPR
          IF (LPRFMML) PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WPR
        ELSEIF (ITYP.EQ.3) THEN
          IF (LERFIIO) ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WPR
          IF (LPRFIIO) PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WPR
        ELSEIF (ITYP.EQ.0) THEN
          IF (LERFPHPHT) 
     .      ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+E0*WPR
          IF (LPRFPHPHT) 
     .      PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WPR
        ENDIF
        IF (MSURFG.GT.0) THEN
          IF (ITYP.EQ.1) THEN
            IF (LERFAAT) ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WPR
            IF (LPRFAAT) PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WPR
          ELSEIF (ITYP.EQ.2) THEN
            IF (LERFMML) ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WPR
            IF (LPRFMML) PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WPR
          ELSEIF (ITYP.EQ.3) THEN
            IF (LERFIIO) ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WPR
            IF (LPRFIIO) PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WPR
          ELSEIF (ITYP.EQ.0) THEN
            IF (LERFPHPHT) 
     .        ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)+E0*WPR
          IF (LPRFPHPHT) 
     .        PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WPR
          ENDIF
        ENDIF
C
C  EITHER: SEMI-TRANSPARENT SURFACE
C
        IF (LTRANS) THEN
          IF (NADSI.GE.1) CALL UPSUSR (WPR,2)
          IF (NADSPC.GE.1) CALL CALC_SPECTRUM (WPR,2)
          RETURN 2
C
C  OR: MIRROR
C
        ELSEIF (ILIIN(MSURF).EQ.3) THEN
          COSI2=-2.*(VELX*CRTX+VELY*CRTY+VELZ*CRTZ)
          VELX=VELX+COSI2*CRTX
          VELY=VELY+COSI2*CRTY
          VELZ=VELZ+COSI2*CRTZ
          IF (NADSI.GE.1) CALL UPSUSR (WPR,2)
          IF (NADSPC.GE.1) CALL CALC_SPECTRUM (WPR,2)
          IF (ICOL.EQ.1) RETURN 3
          RETURN 1
        ENDIF
      ENDIF
C
C   .........................
C   .                       .
C   .  TRANSPARENT SURFACE  .
C   .........................
C
      IF (ILIIN(MSURF).LT.0) THEN
C
C  ONE SIDED FLUX: NEGATIVE COMPONENT
C  IN CASE ILIIN=-3: NET FLUXES HAVE ALREADY BEEN UPDATED ABOVE.
C                    NEED NOT BE UPDATED AGAIN HERE.
C
        IF (SG.GT.0.D0.OR.ILIIN(MSURF).EQ.3) GOTO 90
C
C ITOLD=ITNEW=ITYP
C
        IF (ITYP.EQ.1) THEN
          IF (LERFAAT) ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WPR
          IF (LPRFAAT) PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WPR
        ELSEIF (ITYP.EQ.2) THEN
          IF (LERFMML) ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WPR
          IF (LPRFMML) PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WPR
        ELSEIF (ITYP.EQ.3) THEN
          IF (LERFIIO) ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WPR
          IF (LPRFIIO) PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WPR
        ELSEIF (ITYP.EQ.0) THEN
          IF (LERFPHPHT) 
     .      ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+E0*WPR
          IF (LPRFPHPHT) 
     .      PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WPR
        ENDIF
        IF (MSURFG.GT.0) THEN
          IF (ITYP.EQ.1) THEN
            IF (LERFAAT) ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WPR
            IF (LPRFAAT) PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WPR
          ELSEIF (ITYP.EQ.2) THEN
            IF (LERFMML) ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WPR
            IF (LPRFMML) PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WPR
          ELSEIF (ITYP.EQ.3) THEN
            IF (LERFIIO) ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WPR
            IF (LPRFIIO) PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WPR
          ELSEIF (ITYP.EQ.0) THEN
            IF (LERFPHPHT) 
     .        ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)+E0*WPR
            IF (LPRFPHPHT) 
     .        PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WPR
          ENDIF
        ENDIF
C
90      CONTINUE
        IF (NADSI.GE.1) CALL UPSUSR (WPR,2)
        IF (NADSPC.GE.1) CALL CALC_SPECTRUM (WPR,2)
        RETURN 2
      ENDIF
C
C
100   CONTINUE
C
C   .............................
C   .                           .
C   .  NON TRANSPARENT SURFACE  .
C   .............................
C
      XGENER=0.
C
C  ALL INCIDENT SURFACE TALLIES (INCLUDING SPUTTERING TALLIES, AND
C  INCLUDING CONDITIONAL EXPECTATION CORRECTION, ARE UPDATED.
C  TRANSPARENT, MIRROR AND ABSORBING SURFACES ARE ALSO DONE
C
C  CONDITIONAL EXPECTATION ESTIMATOR: HAS THIS PARTICLE COLLIDED IN THE VOLUME,
C  BEFORE IT HIT THE WALL?
C
      IF (ICOL.EQ.1) RETURN 3
C
C  ..........................................................................
C
C  NOW DEAL WITH REFLECTED AND/OR SPUTTERED PARTICLES.
C  ..........................................................................
C
C  REFLECTION FROM SURFACE
C
C  ......................................................
C  .                                                    .
C  .  REFLECTION MODEL 600--699 FOR INCIDENT MOLECULES  .
C  ......................................................
C
C
600   CONTINUE
C
      IF (ITYP.EQ.2) THEN
C
C  TEST FOR REFLECTION
C
        IF (WEIGHT.GE.WMINS) THEN
C  WITH SUPPRESSION OF ABSORPTION
          WABS=WEIGHT*(1.D0-RECYCT(ISPZ,MSURF))
          IF (WABS.GT.0.D0) THEN
            IF (LSPUMP) SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WABS
            WEIGHT=WEIGHT-WABS
          ENDIF
          IF (WEIGHT.GT.EPS30) GOTO 610
          LGPART=.FALSE.
          RETURN
        ELSEIF (RECYCT(ISPZ,MSURF).NE.1.D0) THEN
C  NO SUPPRESSION OF ABSORPTION
          ZVZ=RANF_EIRENE( )
          IF (ZVZ.LT.RECYCT(ISPZ,MSURF)) GOTO 610
C  ABSORB THIS PARTICLE
          IF (LSPUMP) SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
          LGPART=.FALSE.
          RETURN
        ENDIF
C
610     CONTINUE
C
C  NEW SPECIES: AGAIN MOLECULE
C
C       ITYP=2
        IMOL=ISRT(ISPZ,MSURF)
        IF (IMOL.GT.NMOLI) THEN
          FR2=RANF_EIRENE( )
          DO 621 I=1,NMOLI
            IMOL=I
            IF (FR2.LE.DMOL(IMOL)) GOTO 622
621       CONTINUE
          GOTO 995
622       CONTINUE
        ELSEIF (IMOL.EQ.0) THEN
C  NO THERMAL EMISSION, ABSORB INSTEAD
          IF (LSPUMP) SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
          LGPART=.FALSE.
          RETURN
        ELSEIF (IMOL.LT.0) THEN
C  NOT IN USE          
          GOTO 995
        ENDIF
C
        ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
        E0TERM=EWALL(MSURF)
        IF (E0TERM.GT.0) THEN
C  MONOENERGETIC DISTRIBUTION
          E0=E0TERM
          E0_MEAN=E0
          VEL=RSQDVM(IMOL)*SQRT(E0)
C   AZIMUTAL ANGLE: EQUIDISTRIBUTION
C   POLAR ANGLE: COSINE
          IF (INIV4.EQ.0) CALL FCOSIN
          VX=FC1(INIV4)
          VY=FC2(INIV4)
          VZ=FC3(INIV4)
          INIV4=INIV4-1
          CALL ROTATF (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ)
        ELSEIF (E0TERM.LT.0.D0) THEN
C  SAMPLE FROM MAXWELLIAN FLUX AROUND INNER (!) NORMAL AT TEMP. TW (EV)
          TW=-E0TERM
          CALL VELOCS (TW,0._DP,0._DP,0._DP,0._DP,0._DP,RSQDVM(IMOL),
     .                 CVRSSM(IMOL),
     .                -CRTX,-CRTY,-CRTZ,
     .                 E0,VELX,VELY,VELZ,VEL)
          E0_MEAN=2.*TW
        ELSE
          WRITE (6,*) 'ERROR IN ESCAPE, EXIT CALLED '
          CALL EXIT_OWN(1)
        ENDIF
C
C  ...............................................
C  .                                             .
C  .  REFLECTION MODEL FOR ATOMS OR ATOMIC IONS  .
C  ...............................................
C
      ELSEIF (ITYP.EQ.1.OR.ITYP.EQ.3) THEN
C
C
C  THESE INCIDENT PARTICLES MAY HAVE SPUTTERED AT THIS  SURFACE
C
        SPLFLG=0.
C
        IF (WGHTSP.GT.0..AND.ISSPTP.GT.0) THEN
C  FOLLOW SPUTTERED PARTICLES LATER. PUT THEM INTO STATISTICAL CELLAR
C
          IF (NLEVEL.GE.15) THEN
            WRITE (6,*) 'SPLITTING ABANDONED FOR PART. NO. ',NPANU
            WRITE (6,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL
            GOTO 4712
          ENDIF
          SPLFLG=SPLFLG+1.
          ISPZ=ISSPTP
          ITYP=ISPEZI(ISPZ,-1)
          IPHOT=ISPEZI(ISPZ,0)
          IATM=ISPEZI(ISPZ,1)
          IMOL=ISPEZI(ISPZ,2)
          IION=ISPEZI(ISPZ,3)
          IPLS=ISPEZI(ISPZ,4)
          E0=ESPTP
          WEIGHT=WGHTSP
          VEL=VSPTP
          VELX=VXSPTP
          VELY=VYSPTP
          VELZ=VZSPTP
C  ITOLD.NE.ITNEW=ITYP POSSIBLE
C
          IF (ITYP.EQ.0) THEN
            LOGPHOT(IPHOT,ISTRA)=.TRUE.
            IF (ITOLD.EQ.0) THEN
              IF (LPRFPHPHT)
     .          PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WEIGHT
              IF (LERFPHPHT)
     .          ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFPHPHT)
     .            PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WEIGHT
                IF (LERFPHPHT)
     .            ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)
     .                                   +E0*WEIGHT
              ENDIF
            ENDIF
          ELSEIF (ITYP.EQ.1) THEN
            LOGATM(IATM,ISTRA)=.TRUE.
            IF (ITOLD.EQ.1) THEN
              IF (LPRFAAT) PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WEIGHT
              IF (LERFAAT) 
     .          ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFAAT) 
     .            PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WEIGHT
                IF (LERFAAT)
     .            ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.2) THEN
              IF (LPRFMAT) PRFMAT(IATM,MSURF)=PRFMAT(IATM,MSURF)+WEIGHT
              IF (LERFMAT) 
     .          ERFMAT(IATM,MSURF)=ERFMAT(IATM,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFMAT) 
     .            PRFMAT(IATM,MSURFG)=PRFMAT(IATM,MSURFG)+WEIGHT
              IF (LERFMAT)
     .          ERFMAT(IATM,MSURFG)=ERFMAT(IATM,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.3) THEN
              IF (LPRFIAT) PRFIAT(IATM,MSURF)=PRFIAT(IATM,MSURF)+WEIGHT
              IF (LERFIAT) 
     .          ERFIAT(IATM,MSURF)=ERFIAT(IATM,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFIAT) 
     .            PRFIAT(IATM,MSURFG)=PRFIAT(IATM,MSURFG)+WEIGHT
                IF (LERFIAT)
     .            ERFIAT(IATM,MSURFG)=ERFIAT(IATM,MSURFG)+E0*WEIGHT
              ENDIF
            ENDIF
          ELSEIF (ITYP.EQ.2) THEN
            LOGMOL(IMOL,ISTRA)=.TRUE.
            IF (ITOLD.EQ.1) THEN
              IF (LPRFAML) PRFAML(IMOL,MSURF)=PRFAML(IMOL,MSURF)+WEIGHT
              IF (LERFAML) 
     .          ERFAML(IMOL,MSURF)=ERFAML(IMOL,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFAML) 
     .            PRFAML(IMOL,MSURFG)=PRFAML(IMOL,MSURFG)+WEIGHT
                IF (LERFAML)
     .            ERFAML(IMOL,MSURFG)=ERFAML(IMOL,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.2) THEN
              IF (LPRFMML) PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WEIGHT
              IF (LERFMML) 
     .          ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFMML) 
     .            PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WEIGHT
                IF (LERFMML)
     .            ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.3) THEN
              IF (LPRFIML) PRFIML(IMOL,MSURF)=PRFIML(IMOL,MSURF)+WEIGHT
              IF (LERFIML) 
     .          ERFIML(IMOL,MSURF)=ERFIML(IMOL,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFIML) 
     .            PRFIML(IMOL,MSURFG)=PRFIML(IMOL,MSURFG)+WEIGHT
                IF (LERFIML)
     .            ERFIML(IMOL,MSURFG)=ERFIML(IMOL,MSURFG)+E0*WEIGHT
              ENDIF
            ENDIF
          ELSEIF (ITYP.EQ.3) THEN
            LOGION(IION,ISTRA)=.TRUE.
            IF (ITOLD.EQ.1) THEN
              IF (LPRFAIO) PRFAIO(IION,MSURF)=PRFAIO(IION,MSURF)+WEIGHT
              IF (LERFAIO) 
     .          ERFAIO(IION,MSURF)=ERFAIO(IION,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFAIO) 
     .            PRFAIO(IION,MSURFG)=PRFAIO(IION,MSURFG)+WEIGHT
                IF (LERFAIO)
     .            ERFAIO(IION,MSURFG)=ERFAIO(IION,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.2) THEN
              IF (LPRFMIO) PRFMIO(IION,MSURF)=PRFMIO(IION,MSURF)+WEIGHT
              IF (LERFMIO) 
     .          ERFMIO(IION,MSURF)=ERFMIO(IION,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFMIO) 
     .            PRFMIO(IION,MSURFG)=PRFMIO(IION,MSURFG)+WEIGHT
                IF (LERFMIO)
     .            ERFMIO(IION,MSURFG)=ERFMIO(IION,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.3) THEN
              IF (LPRFIIO) PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WEIGHT
              IF (LERFIIO) 
     .          ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFIIO) 
     .            PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WEIGHT
                IF (LERFIIO)
     .            ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WEIGHT
              ENDIF
            ENDIF
          ENDIF
          IF (NADSI.GE.1) CALL UPSUSR (WEIGHT,2)
          IF (NADSPC.GE.1) CALL CALC_SPECTRUM (WEIGHT,2)
C
C.....................................................................
C  SPLITTING
C
          NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
          DO 533 J=1,NPARTC
            RSPLST(NLEVEL,J)=RPST(J)
533       CONTINUE
          DO 534 J=1,MPARTC
            ISPLST(NLEVEL,J)=IPST(J)
534       CONTINUE
C  NUMBER OF NODES AT THIS LEVEL
          NODES(NLEVEL)=2
C
        ENDIF
C
C  SPLITTING FOR PHYSICAL SPUTTERING DONE

C
        IF (WGHTSC.GT.0..AND.ISSPTC.GT.0) THEN
C
C  CHEMICAL SPUTTERING
C
C  FOLLOW SPUTTERED PARTICLES LATER. PUT THEM INTO STATISTICAL CELLAR
          IF (NLEVEL.GE.15) THEN
            WRITE (6,*) 'SPLITTING ABANDONED FOR PART. NO. ',NPANU
            WRITE (6,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL
            GOTO 4712
          ENDIF
          SPLFLG=SPLFLG+1.
          ISPZ=ISSPTC
          ITYP=ISPEZI(ISPZ,-1)
          IPHOT=ISPEZI(ISPZ,0)
          IATM=ISPEZI(ISPZ,1)
          IMOL=ISPEZI(ISPZ,2)
          IION=ISPEZI(ISPZ,3)
          IPLS=ISPEZI(ISPZ,4)
          E0=ESPTC
          WEIGHT=WGHTSC
          VEL=VSPTC
          VELX=VXSPTC
          VELY=VYSPTC
          VELZ=VZSPTC
C  ITOLD.NE.ITNEW=ITYP POSSIBLE
C
          IF (ITYP.EQ.0) THEN
            LOGPHOT(IPHOT,ISTRA)=.TRUE.
            IF (ITOLD.EQ.0) THEN
              IF (LPRFPHPHT)
     .          PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WEIGHT
              IF (LERFPHPHT)
     .          ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFPHPHT)
     .            PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WEIGHT
                IF (LERFPHPHT)
     .            ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)
     .                                   +E0*WEIGHT
              ENDIF
            ENDIF
          ELSEIF (ITYP.EQ.1) THEN
            LOGATM(IATM,ISTRA)=.TRUE.
            IF (ITOLD.EQ.1) THEN
              IF (LPRFAAT) PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WEIGHT
              IF (LERFAAT) 
     .          ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFAAT) 
     .            PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WEIGHT
                IF (LERFAAT)
     .            ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.2) THEN
              IF (LPRFMAT) PRFMAT(IATM,MSURF)=PRFMAT(IATM,MSURF)+WEIGHT
              IF (LERFMAT) 
     .          ERFMAT(IATM,MSURF)=ERFMAT(IATM,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFMAT) 
     .            PRFMAT(IATM,MSURFG)=PRFMAT(IATM,MSURFG)+WEIGHT
              IF (LERFMAT)
     .          ERFMAT(IATM,MSURFG)=ERFMAT(IATM,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.3) THEN
              IF (LPRFIAT) PRFIAT(IATM,MSURF)=PRFIAT(IATM,MSURF)+WEIGHT
              IF (LERFIAT) 
     .          ERFIAT(IATM,MSURF)=ERFIAT(IATM,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFIAT) 
     .            PRFIAT(IATM,MSURFG)=PRFIAT(IATM,MSURFG)+WEIGHT
                IF (LERFIAT)
     .            ERFIAT(IATM,MSURFG)=ERFIAT(IATM,MSURFG)+E0*WEIGHT
              ENDIF
            ENDIF
          ELSEIF (ITYP.EQ.2) THEN
            LOGMOL(IMOL,ISTRA)=.TRUE.
            IF (ITOLD.EQ.1) THEN
              IF (LPRFAML) PRFAML(IMOL,MSURF)=PRFAML(IMOL,MSURF)+WEIGHT
              IF (LERFAML) 
     .          ERFAML(IMOL,MSURF)=ERFAML(IMOL,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFAML) 
     .            PRFAML(IMOL,MSURFG)=PRFAML(IMOL,MSURFG)+WEIGHT
                IF (LERFAML)
     .            ERFAML(IMOL,MSURFG)=ERFAML(IMOL,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.2) THEN
              IF (LPRFMML) PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WEIGHT
              IF (LERFMML) 
     .          ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFMML) 
     .            PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WEIGHT
                IF (LERFMML)
     .            ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.3) THEN
              IF (LPRFIML) PRFIML(IMOL,MSURF)=PRFIML(IMOL,MSURF)+WEIGHT
              IF (LERFIML) 
     .          ERFIML(IMOL,MSURF)=ERFIML(IMOL,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFIML) 
     .            PRFIML(IMOL,MSURFG)=PRFIML(IMOL,MSURFG)+WEIGHT
                IF (LERFIML)
     .            ERFIML(IMOL,MSURFG)=ERFIML(IMOL,MSURFG)+E0*WEIGHT
              ENDIF
            ENDIF
          ELSEIF (ITYP.EQ.3) THEN
            LOGION(IION,ISTRA)=.TRUE.
            IF (ITOLD.EQ.1) THEN
              IF (LPRFAIO) PRFAIO(IION,MSURF)=PRFAIO(IION,MSURF)+WEIGHT
              IF (LERFAIO) 
     .          ERFAIO(IION,MSURF)=ERFAIO(IION,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFAIO) 
     .            PRFAIO(IION,MSURFG)=PRFAIO(IION,MSURFG)+WEIGHT
                IF (LERFAIO)
     .            ERFAIO(IION,MSURFG)=ERFAIO(IION,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.2) THEN
              IF (LPRFMIO) PRFMIO(IION,MSURF)=PRFMIO(IION,MSURF)+WEIGHT
              IF (LERFMIO) 
     .          ERFMIO(IION,MSURF)=ERFMIO(IION,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFMIO) 
     .            PRFMIO(IION,MSURFG)=PRFMIO(IION,MSURFG)+WEIGHT
                IF (LERFMIO)
     .            ERFMIO(IION,MSURFG)=ERFMIO(IION,MSURFG)+E0*WEIGHT
              ENDIF
            ELSEIF (ITOLD.EQ.3) THEN
              IF (LPRFIIO) PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WEIGHT
              IF (LERFIIO) 
     .          ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WEIGHT
              IF (MSURFG.GT.0) THEN
                IF (LPRFIIO) 
     .            PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WEIGHT
                IF (LERFIIO)
     .            ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WEIGHT
              ENDIF
            ENDIF
          ENDIF
          IF (NADSI.GE.1) CALL UPSUSR (WEIGHT,2)
          IF (NADSPC.GE.1) CALL CALC_SPECTRUM (WEIGHT,2)
C
C.....................................................................
C  SPLITTING
C
          NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
          DO 535 J=1,NPARTC
            RSPLST(NLEVEL,J)=RPST(J)
535       CONTINUE
          DO 536 J=1,MPARTC
            ISPLST(NLEVEL,J)=IPST(J)
536       CONTINUE
C  NUMBER OF NODES AT THIS LEVEL
          NODES(NLEVEL)=2
C
        ENDIF
C
C  SPLITTING FOR CHEMICAL SPUTTERING DONE.
C
C  RESTORE INCIDENT PARTICLE, FOR SURFACE REFLECTION ROUTINE
C
4712    CONTINUE
        IF (SPLFLG.NE.0) THEN
          E0=E0S
          WEIGHT=WEIGHS
          VEL=VELS
          VELX=VELXS
          VELY=VELYS
          VELZ=VELZS
          ISPZ=ISPZS
          LGPART=.FALSE.
        ENDIF
C
        CALL REFLC1 (WMINS,FMASS,FCHAR,NPRT(ISPZ),
     .               ISRF(ISPZ,MSURF),ISRT(ISPZ,MSURF))
        ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
        IF (.NOT.LGPART) THEN
          WEIGHT=0.
        ENDIF
C
C  ............................................
C  .                                          .
C  .  REFLECTION MODEL 700...799 FOR PHOTONS  .
C  ............................................
C
      ELSEIF (ITYP.EQ.0) THEN

        CALL REFLC1_PHOTON (WMINS,FMASS,FCHAR,NPRT(ISPZ),
     .               ISRF(ISPZ,MSURF),ISRT(ISPZ,MSURF))
        ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
        IF (.NOT.LGPART) THEN
          WEIGHT=0.
        ENDIF
      ENDIF
C
C  UPDATE REFLECTED PARTICLE AND ENERGY FLUX
C
      IF (.NOT.LGPART) RETURN
C
C  ITOLD.NE.ITNEW=ITYP POSSIBLE
C
      IF (ITYP.EQ.0) THEN
        LOGPHOT(IPHOT,ISTRA)=.TRUE.
        IF (ITOLD.EQ.0) THEN
          IF (LPRFPHPHT) 
     .      PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WEIGHT
          IF (LERFPHPHT) 
     .      ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFPHPHT) 
     .        PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WEIGHT
            IF (LERFPHPHT) 
     .        ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)+E0*WEIGHT
          ENDIF
        ENDIF
      ELSEIF (ITYP.EQ.1) THEN
        LOGATM(IATM,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          IF (LPRFAAT) PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WEIGHT
          IF (LERFAAT) ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFAAT) PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WEIGHT
            IF (LERFAAT) 
     .        ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.2) THEN
          IF (LPRFMAT) PRFMAT(IATM,MSURF)=PRFMAT(IATM,MSURF)+WEIGHT
          IF (LERFMAT) ERFMAT(IATM,MSURF)=ERFMAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFMAT) PRFMAT(IATM,MSURFG)=PRFMAT(IATM,MSURFG)+WEIGHT
            IF (LERFMAT) 
     .        ERFMAT(IATM,MSURFG)=ERFMAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.3) THEN
          IF (LPRFIAT) PRFIAT(IATM,MSURF)=PRFIAT(IATM,MSURF)+WEIGHT
          IF (LERFIAT) ERFIAT(IATM,MSURF)=ERFIAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFIAT) PRFIAT(IATM,MSURFG)=PRFIAT(IATM,MSURFG)+WEIGHT
            IF (LERFIAT) 
     .        ERFIAT(IATM,MSURFG)=ERFIAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
        ENDIF
      ELSEIF (ITYP.EQ.2) THEN
        LOGMOL(IMOL,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          IF (LPRFAML) PRFAML(IMOL,MSURF)=PRFAML(IMOL,MSURF)+WEIGHT
          IF (LERFAML) ERFAML(IMOL,MSURF)=ERFAML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFAML) PRFAML(IMOL,MSURFG)=PRFAML(IMOL,MSURFG)+WEIGHT
            IF (LERFAML) 
     .        ERFAML(IMOL,MSURFG)=ERFAML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.2) THEN
          IF (LPRFMML) PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WEIGHT
          IF (LERFMML) ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFMML) PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WEIGHT
            IF (LERFMML) 
     .        ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.3) THEN
          IF (LPRFIML) PRFIML(IMOL,MSURF)=PRFIML(IMOL,MSURF)+WEIGHT
          IF (LERFIML) ERFIML(IMOL,MSURF)=ERFIML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFIML) PRFIML(IMOL,MSURFG)=PRFIML(IMOL,MSURFG)+WEIGHT
            IF (LERFIML) 
     .        ERFIML(IMOL,MSURFG)=ERFIML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
        ENDIF
      ELSEIF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          IF (LPRFAIO) PRFAIO(IION,MSURF)=PRFAIO(IION,MSURF)+WEIGHT
          IF (LERFAIO) ERFAIO(IION,MSURF)=ERFAIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFAIO) PRFAIO(IION,MSURFG)=PRFAIO(IION,MSURFG)+WEIGHT
            IF (LERFAIO) 
     .        ERFAIO(IION,MSURFG)=ERFAIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.2) THEN
          IF (LPRFMIO) PRFMIO(IION,MSURF)=PRFMIO(IION,MSURF)+WEIGHT
          IF (LERFMIO) ERFMIO(IION,MSURF)=ERFMIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFMIO) PRFMIO(IION,MSURFG)=PRFMIO(IION,MSURFG)+WEIGHT
            IF (LERFMIO) 
     .        ERFMIO(IION,MSURFG)=ERFMIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
        ELSEIF (ITOLD.EQ.3) THEN
          IF (LPRFIIO) PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WEIGHT
          IF (LERFIIO) ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            IF (LPRFIIO) PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WEIGHT
            IF (LERFIIO) 
     .        ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
        ENDIF
      ENDIF
      IF (NADSI.GE.1) CALL UPSUSR (WEIGHT,2)
      IF (NADSPC.GE.1) CALL CALC_SPECTRUM (WEIGHT,2)
C
      IF  (ITYP.EQ.ITOLD) RETURN 1
      IF ((ITYP.EQ.1.AND.ITOLD.EQ.2).OR.
     .    (ITYP.EQ.2.AND.ITOLD.EQ.1)) RETURN 1
      RETURN
C
995   CONTINUE
      WRITE (6,*) 'SPECIES INDEX OUT OF RANGE IN ESCAPE '
      WRITE (6,*) 'IMOL, MSURF ',IMOL,MSURF
      CALL EXIT_OWN(1)
      END
