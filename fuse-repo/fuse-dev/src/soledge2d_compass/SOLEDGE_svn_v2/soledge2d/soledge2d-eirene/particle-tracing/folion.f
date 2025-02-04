C  MAY05: CALL UPDATE FROM STATIC LOOP WITH IFLAG=4 (RATHER =1)
C         WG. COLL EST. ON 1ST FLIGHT AFTER BIRTH.
C  Sept 05: also vel=velpar before call  to ...col  routines.
!PB 12.01.06: calls to calc_spectrum introduced for cell based spectra
!PB 18.04.06: xstorv=0 in "vacuum region" added
!DR  4.08.06: check v_par=0, otherwise stop trajectory (lable 992)
!DR 10.08.06: cut off Ti with T_vac for collision frequency, for
!             ion tracing in vacuum region
!DR 10.08.06: introduce LCART: TRUE, if velx,vely,velx,vel are cartesian
!                              FALSE,if velx,vely,velx,vel guiding centre
!                                    in this case, cartesian velocity
!                                    is stored in: velxs, velys, velzs, vels
!PB 28.09.06: sg corrected for levgeo=4 and levgeo=5
!DR 09.02.07: not only the direction, but also the magnitute of velocity
!             is reset to full cartesian velocity in subr. NEWFIELD
!PB 22.03.07: LEVGEO=6 --> LEVGEO=10
C
      SUBROUTINE EIRENE_FOLION
C
C     CHARGED PARTICLE, LAUNCHED AT X0,Y0,Z0 IN CELL NRCELL, IPOLG,
C     IPERID, NPCELL, NTCELL, NACELL, NBLOCK, WITH VELOCITY VELX,VELY,VELX
C     IS FOLLOWED.
C     (MODULE: COMPRT.F)
c
c
c
c
c
c
c
C
C  ON INPUT:
C     ITYP=3
C  ON OUTPUT:
C
C           ITYP=0  NEXT GENERATION PHOTON IPHOT IS GENERATED
C           ITYP=1  NEXT GENERATION ATOM IATM IS GENERATED
C           ITYP=2  NEXT GENERATION MOLECULE IMOL IS GENERATED
C           ITYP=4  NO NEXT GENERATION PARTICLE IS GENERATED
C                   (PARTICLE ABSORBED IN BULK ION SPECIES)
C
C  DIFFERENCES FROM SUBR. FOLNEUT:
C    1) MOTION ALONG B (VELPAR,VERPER,....)
C    2) ADDITIONALLY: "FOKKER PLANCK COLLISIONS", ISRFCL=4
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CFPLK
      USE EIRMOD_CLOGAU
      USE EIRMOD_CRAND
      USE EIRMOD_CINIT
      USE EIRMOD_CUPD
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CZT1
      USE EIRMOD_CGEOM
      USE EIRMOD_CTETRA
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CTRIG
      USE EIRMOD_CTRCEI
      IMPLICIT NONE
 
      REAL(DP) :: CFLAG(7,3), DUMT(3), DUMV(3)
      REAL(DP) :: AX(2),v,vv,vx,vy,vz
      REAL(DP) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .            XSTORV2(NSTORV,N2ND+N3RD),
     .            BVEC_1(3), VVEC(3)
      REAL(DP) :: GYRO, COSIN, XLI, YLI, ZLI, FNUI, DIST,
     .          PR, WS, COLTYP, X0ERR, Y0ERR, Z0ERR,
     .          VELXS, VELYS, VELZS, VELS,
     .          PUX, PUY, SG,
     .          VCOS, ZLOG, ZINT1, ZEP1, ZTST, ZINT2,
     .          ZMFP, PN, SH, EIRENE_FPATHI, ZTC,
     .          FNUEQI, XNI, TI,
     .          DELFAC,TIFAC,
     .          SCOS_NEW, XOLD, YOLD
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: ISTS, EIRENE_LEARC2, NCOUS, ICOU, J, JJ, IPL, 
     .           NRCELL_OLD,
     .           ICO, NLI, NLE, NPCELL_OLD, JCOL, NRC, NTCELL_OLD,
     .           NRCOLD, IPLTI, I, IM, IFLAG, ICOUN,NTEST,
     .           EIRENE_LEARC1, IDUM, IFPB
      LOGICAL :: LCNDEXP
C
C  NO CONDITIONAL EXPECTATION ESTIMATORS FOR TEST IONS
C
C  ENERGY LOSS FREQUENCY (LANGER APPROXIMATION) (1/SEC)
C  NUCL.FUS. 22, NO. 6, (1986) P754
      FNUEQI(XNI,TI)=8.8E-8*XNI*TI**(-1.5)
C
C  ALL CELL INDICES MUST BE KNOWN AT THIS POINT
C  TENTATIVELY ASSUME: A NEXT GENERATION PARTICLE WILL BE BORN
C
C  IC_NEUT, IC_ION: COUNTER FOR GENERATIONS WITHIN STATIC LOOP
      IC_ION=IC_NEUT
      LCART=.TRUE.
100   LGPART=.TRUE.
C  FULL CARTESIAN VELOCITY VECTOR VEL,VELX,VELY,VELZ AT THIS POINT
      IF (.NOT.LCART) GOTO 9921
      IC_ION=IC_ION+1
      XGENER=0
      ico=0
      IF (ITYP.EQ.3.AND.(IION.LE.0.OR.IION.GT.NIONI)) GOTO 998
C
C  THE  CELL NUMBER NRCELL, IPOLG, IPERID, NPCELL, NTCELL, NACELL, NBLOCK
C  WAS ALREADY SET IN CALLING SUBROUTINE MCARLO
C
C  IF NLSRFX, SURFACE INDEX MRSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFY, SURFACE INDEX MPSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFZ, SURFACE INDEX MTSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFA, SURFACE INDEX MASURF MUST BE DEFINED AT THIS POINT
C
C  FIND DIRECTION PARALLEL TO B-FIELD
C  I.E. CONVERT CARTESIAN VELOCITY UNIT VECTOR VELX,VELY,VELX INTO
C       PARALLEL AND PERPENDICULAR UNIT VELOCITY COMPONENTES  VELPAR
C
1005  NUPC(1)=NPCELL-1+(NTCELL-1)*NP2T3
      NCELL=NRCELL+NUPC(1)*NR1P2+NBLCKA
      IF (NCELL.GT.NSBOX.OR.NCELL.LT.1) GOTO 991
C  FIND B-FIELD IN CELL NCELL
      CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,0)

C     IF (INDPRO(5) == 8) THEN
C       CALL EIRENE_VECUSR(1,BBX,BBY,BBZ,1)
C     ELSE
C       BBX=BXIN(NCELL)
C       BBY=BYIN(NCELL)
C       BBZ=BZIN(NCELL)
C     END IF
C     BVEC = (/ BBX, BBY, BBZ /)

1003  CONTINUE
      VCOS=VELX*BBX+VELY*BBY+VELZ*BBZ
      IF (ABS(VCOS).LT.EPS30) GOTO 992
      SIGPAR=SIGN(1._DP,VCOS)
      VELXS=VELX
      VELYS=VELY
      VELZS=VELZ
      VELS=VEL
C  USE B-FIELD LINE AS TRAJECTORY
      VLXPAR=SIGPAR*BBX
      VLYPAR=SIGPAR*BBY
      VLZPAR=SIGPAR*BBZ
      VELPAR=ABS(VEL*VCOS)
      VELPER=SQRT(MAX(0._DP,VEL**2 - VELPAR**2))
C  VELOCITY WITH RESPECT TO B-FIELD IS NOW DEFINED:
C  VELPAR: PARALLEL VELOCITY, ABSOLUTE VALUE
C  SIGPAR: SIGN OF PARALLEL VELOCITY WITH RESPECT TO B
C  VELPER: PERPENDICULAR VELOCITY, ALWAYS NON-NEGATIVE
C  VL_PAR: PARALLEL UNIT SPEED VECTOR, VL_PAR=SIG*B
C     VL_PAR=(/VLXPAR,VLYPAR,VLZPAR/)
C
C  SET ION ENERGY = PARALLEL ENERGY OF THE IONIZED TEST PARTICLE
      E0PAR=CVRSSI(IION)*VELPAR*VELPAR
C
1004  CONTINUE
C
C  FOLLOW MOTION OF TEST ION OR "STATIC APPROXIMATION"?
      IF (NFOLI(IION).EQ.-1.AND.IFPATH.EQ.1) GOTO 1001 ! go to static loop
C
C  the particle may be sitting exactly on a surface (nlsrf...=.true.).
C
C  this part is special for ions: due to projection of velocity
C  onto Gyro Center motion (or even onto B-field) the correct
C  angle relative to surface may be lost (e.g. cosin lt 0 may result).
c  Also NINC may be different, depending on whether computed with full
c  or with reduced (guiding centre) velocity
C
C  Hence: Was the correct new cell number NCELL used, in case of nlsrf?
C  Fiddle around a bit with cell number and flight direction in this case.
c  using the reduced (guiding centre) velocity to find orientation
C  relative to surface, and possibly correct side of surface, i.e. cell
c  number
 
 
      IF (NLSRFX) THEN
 
C  PARTICLE IS EXACTLY ON ONE OF THE RADIAL GRID SURFACES (MRSURF)
C  RADIAL CELL NO. NRCELL MAY BE WRONG
C  CHECK ORIENTATION OF PARALLEL MOTION RELATIV TO RADIDAL COORDINATE
C
        NRCELL_OLD=NRCELL
        IF (LEVGEO.EQ.1) THEN
          SG=SIGN(1._DP,VLXPAR)
          IF (SG.LT.0) THEN
            NRCELL=MRSURF-1
          ELSEIF (SG.GT.0) THEN
            NRCELL=MRSURF
          ENDIF
        ELSEIF (LEVGEO.EQ.2) THEN
          PUX= X0-EP1(MRSURF)
          PUY= Y0/ELL(MRSURF)/ELL(MRSURF)
          PN=SQRT(PUX*PUX+PUY*PUY+EPS60)
          PUX=PUX/PN
          PUY=PUY/PN
          SG=VLXPAR*PUX+VLYPAR*PUY
          IF (ABS(SG) .LT. EPS12) THEN
            NLSRFX=.FALSE.
            SH=SIGN(1._DP,SG)*CELDIA(NCELL)*1.D-2
            X0 = X0 + SH*PUX
            Y0 = Y0 + SH*PUY
          END IF
          IF (SG.LT.0) THEN
            NRCELL=NGHPLS(1,MRSURF,NPCELL)
          ELSEIF (SG.GT.0) THEN
            NRCELL=NGHPLS(3,MRSURF,NPCELL)
          ENDIF
        ELSEIF (LEVGEO.EQ.3) THEN
          IFPB = 1
          XOLD = X0
          YOLD = Y0
          IDUM = NPCELL
          SG=VLXPAR*PLNX(MRSURF,NPCELL)+VLYPAR*PLNY(MRSURF,NPCELL)
          DO
            IF (ABS(SG) .LT. EPS12) THEN
              NLSRFX=.FALSE.
              SH=SIGN(1._DP,SG)*CELDIA(NCELL)*1.D-2
              X0 = XOLD + SH*PLNX(MRSURF,NPCELL)*IFPB
              Y0 = YOLD + SH*PLNY(MRSURF,NPCELL)*IFPB
            END IF
            IF (SG.LT.0) THEN
              NRCELL=NGHPLS(1,MRSURF,NPCELL)
            ELSEIF (SG.GT.0) THEN
              NRCELL=NGHPLS(3,MRSURF,NPCELL)
            ELSE
              NRCELL=EIRENE_LEARC1(X0,Y0,Z0,IDUM,MRSURF-1,MRSURF,
     .                             NLSRFX,NLSRFY,NPANU,'FOLION      ')
            ENDIF
            IF (NPCELL == IDUM) EXIT
            IFPB = -1
          END DO
        ELSEIF (LEVGEO.EQ.4) THEN
          SG=VLXPAR*PTRIX(IPOLG,MRSURF)+
     .       VLYPAR*PTRIY(IPOLG,MRSURF)
          IF (ABS(SG) .LT. EPS12) THEN
            SH=SIGN(1._DP,SG)*CELDIA(NCELL)*1.D-8
            X0 = X0  +SH*PTRIX(IPOLG,MRSURF)
            Y0 = Y0  +SH*PTRIY(IPOLG,MRSURF)
            WRITE (iunout,*) 'ON SURFACE IN FOLION, NPANU = ',NPANU
            WRITE (iunout,*) 'and moving parallel to SURFACE'
            WRITE (iunout,*) 'push into suspected cell, sh = ',sh
c dr: I think, if SG gt.0, then NRCELL should be mmodified !!!
            NLSRFX=.FALSE.
          ELSEIF (SG.GT.0) THEN
            NTEST=NCHBAR(IPOLG,MRSURF)
            IF (NTEST.EQ.0) THEN
c  no neighbor. push back into old cell.
              SH=-CELDIA(NCELL)*1.D-8
              WRITE (iunout,*) 'ON SURFACE IN FOLION, NPANU = ',NPANU
              WRITE (iunout,*) 'push back into old cell: sh = ',sh
              WRITE (iunout,*) 'NRCELL = ',NRCELL
              NLSRFX=.FALSE.
c  strictly: particle should be pushed towards COM.
              X0 = X0  +SH*PTRIX(IPOLG,MRSURF)
              Y0 = Y0  +SH*PTRIY(IPOLG,MRSURF)
            ELSE
c  neighbor found. continue in neighbor cell.
              NRCELL=NTEST
              IPOLG=NSEITE(IPOLG,MRSURF)
              MRSURF=NRCELL
            ENDIF
          ELSEIF (SG.LT.0) THEN
C  CONTINUE FLIGHT IN ORIGINAL CELL.
C  NOTHING TO BE DONE
          ENDIF
        ELSEIF (LEVGEO.EQ.5) THEN
          SG=VLXPAR*PTETX(IPOLG,MRSURF)+
     .       VLYPAR*PTETY(IPOLG,MRSURF)+
     .       VLZPAR*PTETZ(IPOLG,MRSURF)
          IF (ABS(SG) .LT. EPS12) THEN
C  TO BE WRITTEN
            WRITE (iunout,*) 'PARALLEL TO SURFACE IN FOLION ',NPANU
            CALL EIRENE_EXIT_OWN(1)
          ELSEIF (SG.GT.0) THEN
            NRCELL=NTBAR(IPOLG,MRSURF)
            IPOLG=NTSEITE(IPOLG,MRSURF)
            MRSURF=NRCELL
          ELSEIF (SG.LT.0) THEN
C  NOTHING TO BE DONE
          ENDIF
        ELSEIF (LEVGEO.EQ.10) THEN
!PB EXPLICITELY ALLOW FOR LEVGEO=10
!PB NOTHING TO BE DONE
        ELSE
          write (iunout,*) 'levgeo in folion  ', levgeo
          write (iunout,*) 'option not ready, exit called'
          call EIRENE_exit_own(1)
        ENDIF
 
        IF (NRCELL.NE.NRCELL_OLD) THEN
          ico=ico+1
          if (ico.le.1) goto 1005
        ENDIF
 
 
      ELSEIF (NLSRFY) THEN
 
 
c  particle is on one of the poloidal grid surfaces (MPSURF)
C  POLOIDAL CELL NO. NPCELL MAY BE WRONG
C  CHECK ORIENTATION OF PARALLEL MOTION RELATIV TO POLOIDAL COORDINATE
C
        NPCELL_OLD=NPCELL
        IF (LEVGEO.EQ.1) THEN
          SG=SIGN(1._DP,VLYPAR)
          IF (SG.LT.0) THEN
            NPCELL=MPSURF-1
          ELSEIF (SG.GT.0) THEN
            NPCELL=MPSURF
          ENDIF
        ELSEIF (LEVGEO.EQ.2.OR.LEVGEO.EQ.3) THEN
          SG=VLXPAR*PPLNX(NRCELL,MPSURF)+VLYPAR*PPLNY(NRCELL,MPSURF)
          IF (SG.LT.0) THEN
            npcell=nghpls(4,nrcell,mpsurf)
            ipolg=npcell
C  ACCOUNT FOR CUTS, PERIODICITY, ETC.
C           mpsurf is correct
          ELSEIF (SG.GT.0) THEN
            npcell=nghpls(2,nrcell,mpsurf)
            ipolg=npcell
C  ACCOUNT FOR CUTS, PERIODICITY, ETC.
            mpsurf=npcell
          ENDIF
        ENDIF
        IF (NPCELL.NE.NPCELL_OLD) THEN
          ico=ico+1
          if (ico.le.1) goto 1005
        ENDIF
 
 
      ELSEIF (NLSRFZ) THEN
 
 
c  particle is on one of the toroidal grid surfaces (MTSURF)
C  TOROIDAL CELL NO. NTCELL MAY BE WRONG
C  CHECK ORIENTATION OF PARALLEL MOTION RELATIV TO POLOIDAL COORDINATE
C
        NTCELL_OLD=NTCELL
C  VLZPAR IS THE RELEVANT VELOCITY COMPONENT, BOTH FOR
C  NLTRZ AND NLTRT OPTION
        SG=SIGN(1._DP,VLZPAR)
        IF (SG.LT.0) THEN
          NTCELL=MTSURF-1
        ELSEIF (SG.GT.0) THEN
          NTCELL=MTSURF
        ENDIF
        IF (NTCELL.NE.NTCELL_OLD) THEN
          ico=ico+1
          if (ico.le.1) goto 1005
        ENDIF
      ENDIF
C
C AT THIS POINT: V_PARALLEL, V_PERP , GYROPHASE, KNOWN
C
      GOTO 1002
C
1001  CONTINUE
      IF (IC_ION.EQ.1.AND.NLTRC.AND.TRCHST)
     .  WRITE (iunout,*) 'TRAJECTORY ENTERS STATIC LOOP, ITYP=', ITYP
 
C***********************************************************************
C  STATIC APPROXIMATION
C  SIMULATE NEXT COLLISION INSTANTANEOUSLY
C***********************************************************************
 
C  WEIGHT TOO SMALL? STOP HISTORY
      IF (WEIGHT.LT.EPS30) THEN
        LGPART=.FALSE.
        RETURN
      ENDIF
C
C  PARTICLE ON SURFACE ?
      IF (NLSRFX.OR.NLSRFY.OR.NLSRFZ.OR.NLSRFA) THEN
C  EMITTED  ?  CALL COLLIDE, AFTER UPDATE
        IF (IC_ION.EQ.1) THEN
C  FIRST ENTRY INTO "STATIC LOOP", ALWAYS: EMITTED FROM FROM SURFACE
C    (CRTXG,....,...): NORMAL RELATIVE TO DEFAULT SETTINGS
C                      NEEDED LATER IF PARTICLE LEAVES STATIC LOOP
C                      VIA STDCOL OR ADDCOL
          CRTXG=CRTX*SCOS
          CRTYG=CRTY*SCOS
          CRTZG=CRTZ*SCOS
          SCOS = SIGN(1.D0,VLXPAR*CRTXG+VLYPAR*CRTYG+VLZPAR*CRTZG)
          SCOS_SAVE = SCOS
          SCOS_NEW  = SCOS
C  INCIDENT DURING STATIC LOOP?  CALL ESCAPE, AFTER UPDATE
        ELSE
          SCOS_NEW = SIGN(1.D0,VLXPAR*CRTXG+VLYPAR*CRTYG+VLZPAR*CRTZG)
        ENDIF
      ELSE
C  PARTICLE NOT ON SURFACE
        SCOS_SAVE = SCOS
        SCOS_NEW  = SCOS
      ENDIF
C
      NCOU=1
      IF (NR1P2 == 0) THEN
        NUPC(1)=0
      ELSE
        NUPC(1)=(NCELL-NRCELL-NBLCKA)/NR1P2
      END IF
C     IF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
        ZMFP=EIRENE_FPATHI(NCELL,CFLAG,1,1)
C     ENDIF
C  XSTOR IN STATIC LOOP:  NOT NEEDED, BECAUSE NCOU=1
C     XSTOR2(:,:,1)=XSTOR(:,:)
C     XSTORV2(:,1)=XSTORV(:)
C  DECIDE TO FOLLOW OR NOT TO FOLLOW THIS TRACK ON BASIS OF MFP
C
C  TO BE WRITTEN
C
      CLPD(1)=ZMFP
      IF (IUPDTE.GE.1) THEN
        IFLAG=4
        CALL EIRENE_UPDION (XSTOR2,XSTORV2,IFLAG)
        CALL EIRENE_CALC_SPECTRUM (WEIGHT,IFLAG,1)
      ENDIF
      ZTC=0.
C  CARRY OUT INELASTIC COLLISION EVENT, DIRECTLY AT PLACE OF BIRTH
      IF (SCOS_SAVE.EQ.SCOS_NEW) THEN
        GOTO 230
      ELSE
C  AT THIS POINT: PARTICLE INCIDENT ON SURFACE, IC_ION GT 1 NECESSARILY
        IF (ILIIN(MSURF).GT.0) THEN
          SCOS=SCOS_NEW
          GOTO 380
        ELSE
          GOTO 230
        END IF
      ENDIF
C
C
1002  CONTINUE
C  NO STATIC APPROXIMATION, FOLLOW MOTION
C
      IF (IC_ION.GT.1.AND.NLTRC.AND.TRCHST)
     .  WRITE (iunout,*) 'TRAJECTORY LEAVES STATIC LOOP, ITYP=',ITYP
      IF (IC_ION.GT.1.AND.
     .   (NLSRFX.OR.NLSRFY.OR.NLSRFZ.OR.NLSRFA)) THEN
C  PARTICLE CONTINUES FROM SURFACE AND FROM PREVIOUS "STATIC LOOP" ?
        IC_ION=0
        IC_NEUT=0
        SCOS_NEW = SIGN(1.D0,VLXPAR*CRTXG+VLYPAR*CRTYG+VLZPAR*CRTZG)
        IF (SCOS_SAVE.NE.SCOS_NEW) THEN
          SCOS=SCOS_NEW
          ZT=0.D0
          TL=0.D0
          IPOLGN=IPOLG
          IF (LCART) THEN
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            VELS=VEL
            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
          ENDIF
          IF (NLSRFA) THEN
            CALL EIRENE_ADDCOL (X0,Y0,Z0,SCOS,*101,*380)
          ELSEIF (NLSRFX) THEN
            IF (LEVGEO.LE.3) THEN
              ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
              MSURFG=NPCELL+(NTCELL-1)*NP2T3
              IF (ILIIN(NLIM+ISTS) .NE. 0)
     .          CALL EIRENE_STDCOL (ISTS,1,SCOS,*101,*380)
            ELSEIF (LEVGEO.EQ.4) THEN
              ISTS=ABS(INMTI(IPOLGN,MRSURF))
              MSURFG=INSPAT(IPOLGN,MRSURF)
              IF (ILIIN(ISTS) .NE. 0)
     .          CALL EIRENE_STDCOL (ISTS,1,SCOS,*101,*380)
            ELSEIF (LEVGEO.EQ.5) THEN
              ISTS=ABS(INMTIT(IPOLGN,MRSURF))
C             MSURFG= ??
              IF (ILIIN(ISTS) .NE. 0)
     .          CALL EIRENE_STDCOL (ISTS,1,SCOS,*101,*380)
            ELSEIF (LEVGEO.EQ.10) THEN
              ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
C             MSURFG= ??
              IF (ILIIN(ISTS) .NE. 0)
     .          CALL EIRENE_STDCOL (ISTS,1,SCOS,*101,*380)
            ENDIF
          ELSEIF (NLSRFY) THEN
            ISTS=INMP2I(IRCELL,MPSURF,ITCELL)
            MSURFG=NRCELL+(NTCELL-1)*NR1P2
            IF (ILIIN(NLIM+ISTS) .NE. 0)
     .        CALL EIRENE_STDCOL (ISTS,2,SCOS,*101,*380)
          ELSEIF (NLSRFZ) THEN
            ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
            MSURFG=NRCELL+(NPCELL-1)*NR1P2
            IF (ILIIN(NLIM+ISTS) .NE. 0)
     .        CALL EIRENE_STDCOL (ISTS,3,SG,*101,*380)
          ENDIF
        ENDIF
      ENDIF
 
C**********************************************************************
C   STATIC LOOP FINISHED. REGULAR PARTICLE TRACKING CONTINUES
C**********************************************************************
 
      IC_ION=0
      IC_NEUT=0
C
C  PARTICLE IN VOLUME OR ON SURFACE BUT NOT FROM "STATIC LOOP"
C
C  EACH TEST ION TRACK STARTS AT THIS POINT, IC_ION=0 HERE
C
101   CONTINUE
C     IF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
C       NLPR=   : NOT AVAILABLE
        NRC=NRCI(IION)
C     ENDIF
C  WEIGHT TOO SMALL? STOP HISTORY
      IF (WEIGHT.LT.EPS30) THEN
        LGPART=.FALSE.
        RETURN
      ENDIF
      ICOL=0
      JCOL=0
C
      ZEP1=RANF_EIRENE( )
      ZLOG=-LOG(ZEP1)
      ZINT1=0.0
      ZINT2=ZINT1
      AX(1)=1.
      AX(2)=1.
      IF (NLTRA) X01=X0+RMTOR
      X00=X0
      Y00=Y0
      Z00=Z0
      Z01=Z0
C  CLEAR WORK VARIABLES AND: CONTINUE FLIGHTS THROUGH TRANSPARENT
C                            SURFACES FROM THIS POINT
104   CONTINUE
      NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
      CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,1)
      NJUMP=0
      DO I=1,NIMINT
        IM=IIMINT(I)
        TIMINT(IM)=0._DP
        IIMINT(I)=0
      END DO
      NIMINT = 0
      TT=1.D30
      TL=1.D30
      TS=1.D30
      ZTST=1.D30                                            
      ZT=0.0
C
      NCOU=1
      NUPC(1)=0
      NCOUNT(1)=1
      NCOUNP(1)=1
      ISRFCL=-1
C
C TL: DISTANCE TO NEXT ADDITIONAL SURFACE
      IF (NCELL.LE.NOPTIM) THEN
        NLI=NLIMII(NCELL)
        NLE=NLIMIE(NCELL)
      ELSE
        NLI=1
        NLE=NLIMI
      ENDIF
      IF (NLI.LE.NLE) THEN
        CALL EIRENE_TIMEA1
     .  (MSURF,NCELL,NLI,NLE,NTCELL,IPERID,X0,Y0,Z0,TIME,
     .               VLXPAR,VLYPAR,VLZPAR,VELPAR,
     .               MASURF,XLI,YLI,ZLI,SG,TL,NLTRC,LCNDEXP)
C       NLPR= :NOT AVAILABLE FOR TEST IONS
        ZTST=TL
        ZDT1=TL
        CLPD(1)=ZDT1
        IF (MASURF.NE.0) ISRFCL=1
      ENDIF
C
C TT: DISTANCE UNTIL NEXT TIMESTEP LIMIT IS REACHED
C     USE VELPAR INSTEAD OF VEL, BECAUSE ORBIT IS COMPUTED WITH VELPAR
C     LATER: VELPAR --> VEL_GC
      IF (LGTIME) THEN
        TT=(DTIMVI-TIME)*VELPAR
        IF (TT.LT.ZTST) THEN
          ZTST=TT
          ZDT1=TT
          CLPD(1)=ZDT1
          ISRFCL=2
        ENDIF
      ENDIF
C
C FNUI: COLLISION FREQUENCY WITH BACKGROUND IONS.
      FNUI=1.D-30
      IF (NRC.GE.0) THEN
        DO IPL=1,NPLSI
          IPLTI=MPLSTI(IPL)
          IF (.NOT.LGVAC(NCELL,IPL))
     .    FNUI=FNUI+FNUEQI(DIIN(IPL,NCELL),TIIN(IPLTI,NCELL))
        ENDDO
      ENDIF
C TAUE: RELAXATION TIME
      TAUE=1./FNUI
C STEPSIZE=0.1*VEL_PARALLEL*TAUE, I.E. 10 COULOMB COLLISIONS PER RELAX.TIME
C TF: DISTANCE UNTIL NEXT COULOMB COLLISION
C DELFAC: INCREASE STEPSIZE AS E0 APPROACHES 1.5 * TI
      TIFAC=MAX(TVAC,TIIN(1,NCELL))
      DELFAC=1.5_DP*TIFAC/ABS(E0-1.5_DP*TIFAC+EPS60)
C  DELTA_T = TAUE*0.1*DELFAC
C  DELTA_S = DELTA_T * VELPAR  ! = TF
C     USE VELPAR INSTEAD OF VEL, BECAUSE ORBIT IS COMPUTED WITH VELPAR
C     LATER: VELPAR --> VEL_GC
      TF=TAUE*VELPAR*0.1*DELFAC
 
      IF (TF.LT.ZTST) THEN
        ZTST=TF
        ZDT1=TF
        CLPD(1)=ZDT1
        ISRFCL=4
      ENDIF
C
C  SCAN OVER RADIAL CELLS
C
210   CONTINUE
C
C  TS:   DISTANCE TO NEXT RADIAL SURFACE OF STANDARD MESH
C  ZDT1: DISTANCE TRAVELLED IN CURRENT RADIAL CELL
C
C  USE PARALLEL VELOCITY, I.E., COMPUTE PARALLEL DISTANCES IN GRID
C  THUS ZT,TS,ZTST,ZDT1,CLPD ETC. ARE PARALLEL DISTANCES
C
      IF (ITIME.EQ.1) THEN
        IF (LCART) THEN
          VELXS=VELX
          VELYS=VELY
          VELZS=VELZ
          VELS =VEL
          VELX=VLXPAR
          VELY=VLYPAR
          VELZ=VLZPAR
          VEL =VELPAR
          LCART=.FALSE.
        ENDIF
 
        IF (NLRAD) THEN
          CALL EIRENE_TIMER(TS)
C
          IF (TL.LT.TS.OR.TT.LT.TS.OR.TF.LT.TS) THEN
            MRSURF=0
            IPOLGN=0
C  COLLISION WITH ADDITIONAL SURFACE
            IF (TL.LE.TT.AND.TL.LE.TF) THEN
              ZDT1=TL-ZT
              TL=ZT+ZDT1
              ZTST=TL
              ISRFCL=1
C  COLLISION WITH TIME SURFACE
            ELSEIF (TT.LT.TL.AND.TL.LE.TF) THEN
              ZDT1=TT-ZT
              TT=ZT+ZDT1
              ZTST=TT
              ISRFCL=2
C  FOKKER PLANCK COLLISION
            ELSEIF (TF.LT.TL.AND.TF.LE.TT) THEN
              ZDT1=TF-ZT
              TF=ZT+ZDT1
              ZTST=TF
              ISRFCL=4
            ENDIF
          ELSE
C  COLLISION WITH RADIAL SURFACE
            ISRFCL=0
            ZDT1=TS-ZT
            ZTST=TS
          ENDIF
        ENDIF
C
        NCOU=1
        NUPC(1)=0
        CLPD(1)=ZDT1
        NCOUNT(1)=1
        NCOUNP(1)=1
C
        IF (NLTOR.OR.NLTRA) THEN
          CALL EIRENE_TIMET (ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C
        IF (NLPOL) THEN
          CALL EIRENE_TIMEP(ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C
        IF (ZDT1.LE.0.D0) GOTO 990
        IF (.NOT.LCART) THEN
          VELX=VELXS
          VELY=VELYS
          VELZ=VELZS
          VEL =VELS
          LCART=.TRUE.
        ENDIF
 
      ENDIF
      IF (ZTST.GE.1.D30) GOTO 990
C
C  LOCAL MEAN FREE PATH
C  USE PARALLEL VELOCITY, I.E., COMPUTE PARALLEL MFP
C  BECAUSE CLPD IS THE PARALLEL DISTANCE IN EACH CELL (EXCLUD. GYRO)
C  ETC.. E.G LAMBDA(PARALLEL)=VEL(PARALLEL)/SIGV.
C  THE COLLISION FREQUENCY SIGV, HOWEVER, MUST BE COMPUTED USING THE
C  FULL TEST ION VELOCITY VECTOR, BECAUSE IT MAY DEPEND UPON THE RELATIV
C  INTERACTION ENERGY: TO BE WRITTEN
C  FOR INTERACTIONS WITH ELECTRONS THIS IS IRRELEVANT
C
      IF (IFPATH.NE.1.OR.NRC.LT.0) THEN
        XSTORV(:)=0.D0
        DO 214 J=1,NCOU
          JJ=J
          XSTOR2(:,:,J)=0.D0
          XSTORV2(:,J)=0.D0
          ZMFP=1.D10
          IF (NLPOL) NPCELL=NCOUNP(J)
          IF (NLTOR) NTCELL=NCOUNT(J)
C         VEL=VELS
          GOTO 213
214     CONTINUE
      ELSE
        IF (LCART) THEN
          VELXS=VELX
          VELYS=VELY
          VELZS=VELZ
          VELS =VEL
          VELX=VLXPAR
          VELY=VLYPAR
          VELZ=VLZPAR
          VEL =VELPAR
          LCART=.FALSE.
        ENDIF
        DO 212 J=1,NCOU
          JJ=J
          NCELL=NRCELL+NUPC(J)*NR1P2+NBLCKA
          ZMFP=EIRENE_FPATHI(NCELL,CFLAG,J,NCOU)
          IF (NCOU.GT.1) THEN
            XSTOR2(:,:,J)=XSTOR(:,:)
            XSTORV2(:,J)=XSTORV(:)
          ENDIF
C  UPDATE INTEGRAL
          ZINT1=ZINT1+CLPD(J)*ZMFPI
C         IF (.NOT.NLPR) THEN
CCC         IF (ZINT1.GE.ZLOG) THEN
              IF (NLPOL) NPCELL=NCOUNP(J)
              IF (NLTOR) NTCELL=NCOUNT(J)
              VELX=VELXS
              VELY=VELYS
              VELZ=VELZS
              VEL =VELS
              LCART=.TRUE.
              GOTO 213
CCC         ENDIF
C  THESE NEXT TWO LINES CAN NEVER BE REACHED, BECAUSE ONLY ONE
C  CELL FOR EACH TRACK OF IONS (DISTINCT FROM FOLNEUT).
C  THEN (AT THE LATEST): ROTATION OF VELOCITY DUE TO NEW B-FIELD
C           ZINT2=ZINT1
C           ZT=ZT+CLPD(J)
C         ELSEIF (JCOL.EQ.0) THEN
C   CONDITIONAL EXPECTATION ESTIMATOR FOR TEST IONS: TO BE WRITTEN
C         ENDIF
212     CONTINUE
        VELX=VELXS
        VELY=VELYS
        VELZ=VELZS
        VEL =VELS
        LCART=.TRUE.
      ENDIF
C
213   CONTINUE
      NCOUS=NCOU
      NCOU=JJ
 
CCC  IF NO COLLISION, THEN: ENFORCE ONLY ONE STEP AT A TIME
      IF (ZINT1.LT.ZLOG.AND.NCOUS.GT.1) THEN
        MRSURF=0
        MPSURF=0
        MTSURF=0
        MASURF=0
        ISRFCL=0
        NINCX=0
        NINCY=0
        NINCZ=0
      ENDIF
CCC
C
C  CHECK FOR EVENT
C
C     IF (NLPR)    ......
      IF (ZINT1.GE.ZLOG) GO TO 220
C
      ZINT2=ZINT1
      ZT=ZTST
C
C  RESET CLPD TO REAL PATH LENGTH OF GYRO MOTION
      DO 217 ICOU=1,NCOU
        CLPD(ICOU)=CLPD(ICOU)*VEL/VELPAR
217   CONTINUE
C
C  UPDATE CONTRIBUTION TO VOLUME AVERAGED ESTIMATORS
C
      IF (IUPDTE.GE.1) THEN
        CALL EIRENE_UPDION(XSTOR2,XSTORV2,3)
        CALL EIRENE_CALC_SPECTRUM (WEIGHT,3,1)
      ENDIF
C
C  STOP TRACK ?
C
CDR: ALLE DISTANZEN IN ...COL routines sind parallele distanzen
CDR: Daher auch wg. x = x + dist/vel  parallele geschwindigkeiten.
      IF (LCART) THEN
        VELXS=VELX
        VELYS=VELY
        VELZS=VELZ
        VELS =VEL
        VELX=VLXPAR
        VELY=VLYPAR
        VELZ=VLZPAR
        VEL =VELPAR
        LCART=.FALSE.
      ENDIF
      IF (ISRFCL.EQ.1) CALL EIRENE_ADDCOL(XLI,YLI,ZLI,SG,*104,*380)
      IF (ISRFCL.EQ.2) CALL EIRENE_TIMCOL(AX(2),         *104,*800)
      IF (ISRFCL.EQ.3) CALL EIRENE_TORCOL(               *104)
      IF (ISRFCL.EQ.4) CALL EIRENE_FPKCOL(               *104,*100)
      VELX=VELXS
      VELY=VELYS
      VELZ=VELZS
      VEL =VELS
      LCART=.TRUE.
C
C  NO, CONTINUE TRACK
C
216   CONTINUE
C
C  NEXT CELL - CHECK FOR ESCAPE OR NON DEFAULT ACTING STANDARD SURFACE
C
      IF (LEVGEO.LE.3) THEN
C
        ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
        IF (NLRAD.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCX)
          NLSRFX=.TRUE.
          MSURFG=NPCELL+(NTCELL-1)*NP2T3
          IF (LCART) THEN
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            VELS =VEL
            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
          ENDIF
          IF (ILIIN(NLIM+ISTS) .NE. 0) CALL EIRENE_STDCOL
     .  (ISTS,1,SG,*104,*380)
        ENDIF
        ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
        IF (NLTOR.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCZ)
          NLSRFZ=.TRUE.
          MSURFG=NRCELL+(NPCELL-1)*NR1P2
          IF (LCART) THEN
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            VELS =VEL
            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
          ENDIF
          IF (ILIIN(NLIM+ISTS) .NE. 0) CALL EIRENE_STDCOL
     .  (ISTS,3,SG,*104,*380)
        ENDIF
        ISTS=INMP2I(IRCELL,MPSURF,ITCELL)
        IF (NLPOL.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCY)
          NLSRFY=.TRUE.
          MSURFG=NRCELL+(NTCELL-1)*NR1P2
          IF (LCART) THEN
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            VELS =VEL
            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
          ENDIF
          IF (ILIIN(NLIM+ISTS) .NE. 0) CALL EIRENE_STDCOL
     .  (ISTS,2,SG,*104,*380)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.4) THEN
        ISTS=ABS(INMTI(IPOLGN,MRSURF))
        IF (NLRAD.AND.ISTS.NE.0) THEN
!pb          SG=ISIGN(1,NINCX)
          NLSRFX=.TRUE.
          MSURFG=INSPAT(IPOLGN,MRSURF)
          IF (LCART) THEN
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            VELS =VEL
            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
          ENDIF
          SG=SIGN(1._DP,VELX*PTRIX(IPOLGN,MRSURF)+
     .                  VELY*PTRIY(IPOLGN,MRSURF))
          IF (ILIIN(ISTS) .NE. 0) CALL EIRENE_STDCOL
     .  (ISTS,1,SG,*104,*380)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.5) THEN
        ISTS=ABS(INMTIT(IPOLGN,MRSURF))
        IF (NLRAD.AND.ISTS.NE.0) THEN
          NLSRFX=.TRUE.
          IF (LCART) THEN
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            VELS =VEL
            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
          ENDIF
          SG=SIGN(1._DP,VELX*PTETX(IPOLGN,MRSURF)+
     .                  VELY*PTETY(IPOLGN,MRSURF)+
     .                  VELZ*PTETZ(IPOLGN,MRSURF))
          IF (ILIIN(ISTS) .NE. 0) CALL EIRENE_STDCOL
     .  (ISTS,1,SG,*104,*380)
        ENDIF
C
      ELSEIF (LEVGEO.EQ.10) THEN
        ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
        IF (NLRAD.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCX)
          NLSRFX=.TRUE.
          IF (LCART) THEN
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            VELS =VEL
            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
          ENDIF
          IF (ILIIN(ISTS) .NE. 0) CALL EIRENE_STDCOL
     .  (ISTS,1,SG,*104,*380)
        ENDIF
      ENDIF
C
C
      NRCELL=NRCELL+NINCX
      IF (NRCELL.GT.NR1STM.OR.NRCELL.LT.1) GOTO 990
C
CDR: SPLITTING AND COND.EXP.EST. NOT AVAILABLE FOR TEST IONS
C
C  CHECK IF WE HAVE ENCOUNTERED A SPLITTING ZONE
C     IF (NLSPLT(MRSURF).AND.NLEVEL.LT.MAXLEV.AND.ICOL.EQ.0) GOTO 330
C
C  SWITCH OFF CONDITIONAL EXP. ESTIMATOR ?
C     IF (AX(2).LT.WMINC) THEN
C       IF (ICOL.EQ.1) GOTO 512
C  NO COLLISION YET; RESTART AGAIN WITH COND. EXP. ESTIMATOR
C                    IN NEW CELL
C       AX(1)=1.
C       AX(2)=1.
C       JCOL=0
C     ENDIF
CCC
      ZTC=CLPD(1)*VELPAR/VEL
      GOTO 2211
CCC
CCC   GOTO 210
C
C  POINT OF COLLISION  220---240
C
220   CONTINUE
C
      CLPD(NCOU)=(ZLOG-ZINT2)*ZMFP
      ZTC=ZT+CLPD(NCOU)
C  RESET CLPD TO REAL PATH LENGTH OF GYRO MOTION
      DO 221 ICOU=1,NCOU
        CLPD(ICOU)=CLPD(ICOU)*VEL/VELPAR
221   CONTINUE
      IF (IUPDTE.GE.1) THEN
        CALL EIRENE_UPDION (XSTOR2,XSTORV2,4)
        CALL EIRENE_CALC_SPECTRUM (WEIGHT,4,1)
      ENDIF
2211  continue
      X0=X0+VLXPAR*ZTC
      Y0=Y0+VLYPAR*ZTC
      Z0=Z0+VLZPAR*ZTC
      TIME=TIME+ZTC/VELPAR
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'FOLION 2     ')
      ELSEIF (NLFEM) THEN
        IPOLG=0
      ELSEIF (NLTET) THEN
        IPOLG=0
      ENDIF
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      NLSRFA=.FALSE.
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      MSURF=0
      IF (NLTRA) PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
C
CCC

C  DELTA EVENT AT CELL BOUNDARY

      IF (ZINT1.LT.ZLOG) THEN
        IF (NINCX.NE.0) THEN
          NLSRFX=.TRUE.
          IF (LEVGEO < 4) THEN
            MRSURF=NRCELL
            IF (NINCX.EQ.-1) MRSURF=NRCELL+1
          ELSEIF (LEVGEO == 4) THEN
            NRCOLD=NRCELL-NINCX
            MRSURF=NCHBAR(IPOLGN,NRCOLD)
            IPOLG=NSEITE(IPOLGN,NRCOLD)
          ELSEIF (LEVGEO == 5) THEN
            NRCOLD=NRCELL-NINCX
            MRSURF=NTBAR(IPOLGN,NRCOLD)
            IPOLG=NTSEITE(IPOLGN,NRCOLD)
          ELSEIF (LEVGEO == 10) THEN
!PB EXPLICITELY ALLOW FOR LEVGEO=10
!PB NOTHING DONE FOR DELTA EVENT AT CELL BOUNDARY
          ELSE
            WRITE (iunout,*) 'DELTA EVENT AT CELL BOUNDARY NOT READY '
            WRITE (iunout,*) 'FOR LEVGEO=10 IN SUBR. FOLION. '
            CALL EIRENE_EXIT_OWN(1)
          END IF
        ELSEIF (NINCZ.NE.0) THEN
          NLSRFZ=.TRUE.
          NTCELL=KUPC(1)+NINCZ
          IF (NINCZ == 1) THEN
            MTSURF=NTCELL
          ELSEIF (NINCZ.EQ.-1) THEN
            MTSURF=NTCELL+1
          ENDIF
        ELSEIF (NINCY.NE.0) THEN
          NLSRFY=.TRUE.
          IF (LEVGEO.EQ.1) THEN
            NPCELL=JUPC(1)+NINCY
            IF (NINCY == 1) THEN
              MPSURF=NPCELL
            ELSEIF (NINCY.EQ.-1) THEN
              MPSURF=NPCELL+1
            ENDIF
          ELSEIF (LEVGEO.LE.3) THEN
            MPSURF=LUPC(1)
            IF (MUPC(1).EQ.1) NPCELL=NGHPLS(2,NRCELL,MPSURF)
            IF (MUPC(1).NE.1) NPCELL=NGHPLS(4,NRCELL,MPSURF)
C  PERIODICITY FOR LEVGEO=2 (TO BE WRITTEN INTO MORE GENERAL TERMS)
            IF (NPCELL.EQ.0.AND.LEVGEO.EQ.2) THEN
              WRITE (iunout,*) 'should not be here '
              MPSURF=NP2ND
              NPCELL=NP2NDM
            ELSEIF (NPCELL.EQ.NP2ND.AND.LEVGEO.EQ.2) THEN
              WRITE (iunout,*) 'should not be here '
              MPSURF=1
              NPCELL=1
            ENDIF
            IF (LEVGEO.LE.3.AND.NLPOL) THEN
              IPOLG=NPCELL
            ELSEIF (LEVGEO.EQ.3.AND..NOT.NLPOL) THEN
              IPOLG=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'FOLION neu   ')
            ENDIF
          ENDIF
        ELSE
          IF (LEVGEO /= 10) GOTO 994
        ENDIF
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,19)
        NUPC(1)=NPCELL-1+(NTCELL-1)*NP2T3
        NCELL=NRCELL+NUPC(1)*NR1P2+NBLCKA
C  DELTA COLLISION AT SURFACE DONE, NEW CELL FOUND

C  FIND NEW B-FIELD, NEW REDUCED (GC) VELOCITY
        CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,1)

C  FIND NEW MAGNETIC FIELD
c       IF (INDPRO(5) == 8) THEN
c         CALL VECUSR(1,BBX,BBY,BBZ,1)
c       ELSE
c         BBX=BXIN(NCELL)
c         BBY=BYIN(NCELL)
c         BBZ=BZIN(NCELL)
c       END IF
c       BVEC = (/ BBX, BBY, BBZ /)
C  RETAIN VEL, V_PARALLEL, V_PERP, SIGPAR,
C  SAMPLE PHASE, AND FIND NEW CARTESIAN VX,VY,VZ (SAME VEL=VELS)
c       VLXPAR=SIGPAR*BBX
c       VLYPAR=SIGPAR*BBY
c       VLZPAR=SIGPAR*BBZ
c       VELPAR=ABS(VELS*VCOS)*(1.D0-EPS12)
c       E0PAR=CVRSSI(IION)*VELPAR*VELPAR
c       VELPER = SQRT(VELS**2 - VELPAR**2)
C  NEW GYRO PHASE
c       GYRO=RANF_EIRENE()*PI2A
C  BACK TO CARTESIAN COORDIANTES
c       CALL B_PROJI (BVEC,BVEC_1,VVEC,SIGPAR*VELPAR,VELPER,GYRO)
c       VELX = VVEC(1)
c       VELY = VVEC(2)
c       VELZ = VVEC(3)
c       VEL=VELS
c       LCART=.TRUE.

        ICO = 0
        GOTO 1004
      ENDIF
CCC
C
230   CONTINUE
C
C  PRE COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL EIRENE_UPCUSR(WS,1)
        CALL EIRENE_CALC_SPECTRUM (WS,1,1)
      ENDIF
C
C
C  TEST FOR CORRECT CELL NUMBER AT COLLISION POINT
C  KILL PARTICLE, IF TOO LARGE ROUND OFF ERRORS DURING
C  PARTICLE TRACING
C
      IF (NLTEST) CALL EIRENE_CLLTST(*997)
C
C  SAMPLE FROM COLLISION KERNEL FOR TEST IONS
C  FIND NEW WEIGHT, SPECIES INDEX, VELOCITY AND RETURN
C
      CALL EIRENE_COLION(CFLAG,COLTYP,DIST)
      ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C  POST COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL EIRENE_UPCUSR(WS,2)
        CALL EIRENE_CALC_SPECTRUM (WS,2,1)
      ENDIF
C
      IF (COLTYP.EQ.2.) GOTO 700
C
      GOTO 100
C
C  SIMULATION OF COLLISION EVENT FINISHED
C
C
C
C   INCIDENT ONTO SURFACE
C
380   CONTINUE
C
C   REFLECTION FROM  SURFACE
C   USE FULL VELOCITY, NOT ONLY THE PARALLEL VELOCITY.
C   THIS IS DONE BY SAMPLING THE GYRO-PHASE IN SUBR. NEWFIELD
C   REJECT THOSE GYROPHASES WHICH WOULD LEAD TO NEGATIVE ANGLE OF INCIDENCE
C
      IF (.NOT.LCART) THEN
        NUPC(1)=NPCELL-1+(NTCELL-1)*NP2T3
        NCELL=NRCELL+NUPC(1)*NR1P2+NBLCKA
        ICOUN=0
        DO
cdr     write (iunout,*) 'before newfield ', vel,velx,vely,velz,e0
cdr     write (iunout,*) 'before newf., save ', vels,velxs,velys,velzs
cdr     write (iunout,*) 'par,per ',velpar,velper,
cdr  .      sqrt(velpar*velpar+velper*velper)
          CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,2)
          COSIN=VELX*CRTX+VELY*CRTY+VELZ*CRTZ
C  DOES THE PARTICLE SPEED POINT TOWARDS THE SURFACE
          IF (COSIN.GT.0.) EXIT
C  NO, TRY NEXT GYRO PHASE
          ICOUN=ICOUN+1
          IF (ICOUN.EQ.100) THEN
            WRITE (IUNOUT,*) 'PARTICLE KILLED AT SURFACE IN FOLION'
            WRITE (IUNOUT,*) 'NPANU ',NPANU
            WRITE (IUNOUT,*) 'VELPER,VELPAR ',VELPER,VELPAR
            LGPART=.FALSE.
            WEIGHT=0.
            RETURN
          ENDIF
 
        ENDDO
C  NOW A PARTICLE WITH FULL CARTESIAN VELOCITY VECTOR (LCART=0)
C  IS SET. ITS SPEED VECTOR POINTS TOWARDS THE SURFACE (COSIN.GT.0)
      ENDIF
C
C  UPDATE EFFLUXES ONTO SURFACE AND REFLECT PARTICLE
      PR=1.
      IF (ILIIN(MSURF).LE.-2) PR=SG
C
C  FOR NONTRANSPARENT SURFACES:
C  ACCELERATION IN SHEATH IS DONE IN SUBR. ESCAPE
C
      CALL EIRENE_ESCAPE(PR,SG,*100,*104,*996)
      RETURN
C
C   100: START NEW ION TRACK
C   104: CONTINUE THIS TRACK, TRANSPARENT SURFACE IS CROSSED
C
C
700   CONTINUE
C  REGULAR STOP IN SUBR. FOLION, CONTINUE IN SUBR. MCARLO
      RETURN
C
800   CONTINUE
C  REGULAR STOP IN SUBR. FOLION, STOP HISTORY, CENSUS ARRAY FULL
C     IF (ICOL.EQ.1.AND..NOT.LGLAST) GOTO 512
      LGPART=.FALSE.
      WEIGHT=0.
      RETURN
C
990   CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('ERROR IN FOLION,  ZDT1 OR NRCELL OUT OF RANGE  ')
      CALL EIRENE_MASAGE
     .  ('PARTICLE IS KILLED                            ')
      WRITE (iunout,*) 'NPANU,NRCELL,ZDT1,ZTST ',NPANU,NRCELL,ZDT1,ZTST
      WRITE (iunout,*) 'TL,TS,ZINT1,ZLOG ',TL,TS,ZINT1,ZLOG
      GOTO 995
991   CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('ERROR IN FOLION,  NCELL OUT OF RANGE            ')
      CALL EIRENE_MASAGE
     .  ('PARTICLE IS KILLED                            ')
      WRITE (iunout,*) 'NPANU,NCELL,NRCELL,NPCELL,NTCELL '
      WRITE (iunout,*)  NPANU,NCELL,NRCELL,NPCELL,NTCELL
      GOTO 995
992   CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('ERROR IN FOLION,  PROJECTION TO V_PAR, V_PERP   ')
      CALL EIRENE_MASAGE
     .  ('PROBABLY ILL DEFINED B-FIELD WRT. PARTICLE SPEED')
      WRITE (iunout,*) 'BBX,BBY,BBZ ',BBX,BBY,BBZ
      ZT=0.
      GOTO 9951
9921  CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('ERROR IN FOLION,  LCART HAS WRONG VALUE       ')
      WRITE (iunout,*) 'npanu,lcart ',npanu,lcart
      IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
      GOTO 999
993   CALL EIRENE_MASAGE
     .  ('ERROR IN FOLION,  NO PARTICLE TRACING BUT     ')
      CALL EIRENE_MASAGE
     .  ('IFPATH.NE.1. PARTICLE IS KILLED               ')
      WRITE (iunout,*) 'IION ',IION
      GOTO 999
C
994   CALL EIRENE_MASAGE
     .  ('ERROR IN FOLION,  AT SURFACE DELTA EVENT      ')
      WRITE (iunout,*) 'IION,NPANU ',IION,NPANU
      GOTO 999
C
995   WRITE (iunout,*) 'MRSURF,MPSURF,MTSURF,MASURF ',
     .                  MRSURF,MPSURF,MTSURF,MASURF
9951  X0ERR=X0+ZT*VELX
      Y0ERR=Y0+ZT*VELY
      Z0ERR=Z0+ZT*VELZ
      IF (NLTRC) THEN
        CALL EIRENE_CHCTRC(X0ERR,Y0ERR,Z0ERR,16,18)
      ELSE
        WRITE (iunout,*) 'X0,Y0,Z0,ZT ',X0,Y0,Z0,ZT
        WRITE (iunout,*) 'VELX,VELY,VELZ ',VELX,VELY,VELZ
        WRITE (iunout,*) 'X0ERR,Y0ERR,Z0ERR ',X0ERR,Y0ERR,Z0ERR
      ENDIF
      GOTO 999
996   CALL EIRENE_MASAGE
     .  ('ERROR IN FOLION, COND. EXP. ESTIM. NOT IN USE ')
      GOTO 999
997   CALL EIRENE_MASAGE
     .  ('ERROR IN FOLION,  DETECTED IN SUBR. CLLTST    ')
      CALL EIRENE_MASAGE
     .  ('PARTICLE IS KILLED                            ')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST
      IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
      GOTO 999
998   WRITE (iunout,*) 'ERROR IN FOLION, SPECIES INDEX OUT OF RANGE '
      WRITE (iunout,*) ' NPANU,IION ',NPANU,IION
      GOTO 999
C
999   PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      WEIGHT=0.
      CALL EIRENE_LEER(1)
      RETURN
      END
 
      SUBROUTINE EIRENE_NEWFIELD(X,Y,Z,VELS,IND)
C  FIND NEW MAGNETIC FIELD AT NEW POINT X,Y,Z IN CELL NCELL
C
C  IF (IND.GE.1) ALSO PROVIDE REDUCED (GC) VELOCITY
C    RETAIN V_PARALLEL, V_PERP
C    CHECKS DONE THAT VELPER AND VERPAR ARE PRESERVED, CHECKS REMOVED.

C  IF (IND.GE.2) ALSO PROVIDE NEW CARTESIAN VELOCITY
C  BY SAMPLING THE GYRO PHASE, AND A COORDINATE TRANSFORMATION IN
C  VEL-SPACE.
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CFPLK
      USE EIRMOD_COMPRT
      USE EIRMOD_CRAND
      USE EIRMOD_CINIT
      IMPLICIT NONE
      REAL(DP), EXTERNAL :: RANF_EIRENE
      REAL(DP), INTENT(IN) :: X,Y,Z,VELS
      REAL(DP) :: BVEC_1(3), VVEC(3), GYRO
      INTEGER :: IND
 
      IF (INDPRO(5) == 8) THEN
        CALL EIRENE_VECUSR(1,BBX,BBY,BBZ,1)
      ELSE
        BBX=BXIN(NCELL)
        BBY=BYIN(NCELL)
        BBZ=BZIN(NCELL)
      END IF
      BVEC = (/ BBX, BBY, BBZ /)

      IF (IND.LT.1) RETURN

C  RETAIN VEL, V_PARALLEL, V_PERP, SIGPAR,
C  SAMPLE PHASE, AND FIND NEW CARTESIAN VELX,VELY,VELZ (SAME VEL=VELS)
      VLXPAR=SIGPAR*BBX
      VLYPAR=SIGPAR*BBY
      VLZPAR=SIGPAR*BBZ
      VELX = VLXPAR
      VELY = VLYPAR
      VELZ = VLZPAR
      VEL  = VELPAR
      LCART=.FALSE.

      IF (IND.LT.2) RETURN
                                            
C  NEW GYRO PHASE
      GYRO=RANF_EIRENE()*PI2A
C  BACK TO CARTESIAN COORDIANTES
      CALL EIRENE_B_PROJI (BVEC,BVEC_1,VVEC,SIGPAR*VELPAR,VELPER,GYRO)
      VELX = VVEC(1)
      VELY = VVEC(2)
      VELZ = VVEC(3)
      VEL  = VELS
      LCART=.TRUE.
      RETURN
      END
