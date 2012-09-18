c     Krieger IPP 2012 - fixed INTEL name conflict: ranf->ranf_eirene
      FUNCTION EMAXW(TI,XMPER,YMPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  MEAN ENERGY OF AN ION, GIVEN A TRUNCATED MAXWELLIAN FLUX DENSITY
C  AT TI (EV) SHIFTED
C  BY A DRIFT WITH NORMALIZED VELOCITY XMPER PERPENDICULAR
C  ONTO AND YMPAR PARALLEL TO THE TARGET SURFACE
C  IF VDRIFT IS A VELOCITY, THEN XMPER (YMPAR) IS
C  VDRIFT/VTERM WITH VTERM=SQRT(2*TI/RMASS)
C  IN EIRENE UNITS: VDRIFT: CM/SEC
C                   VTERM : SQRT(2*TI/RMASS)*CVEL2A (CM/SEC)
C                           (WITH TI (EV), RMASS (AMU) CVEL2A=9.8226 E5)
      INCLUDE 'PARMMOD'
      INCLUDE 'CCONA'
      IF (XMPER.NE.0.D0) THEN
        XMM=XMPER
        XMM2=XMPER*XMPER
        YMM2=YMPAR*YMPAR
        ER1=EXP(-XMM2)/SQRT(PIA)
        ER2=(1.+ERF(XMM))*XMM
        FACTOR=((XMM2+2.+YMM2)*ER1+(XMM2+2.5+YMM2)*ER2)/(ER1+ER2)
      ELSE
        YMM2=YMPAR*YMPAR
        FACTOR=2.+YMM2
      ENDIF
      EMAXW=FACTOR*TI
      RETURN
      END
C
      FUNCTION SHEATH(TE,DP,VP,NZP,GAMMA,CUR,NP,MS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  JAN.93: NEW VERSION, NOW CORRECTLY ACCOUNTING
C          FOR MULTISPECIES PLASMA BACKGROUND
C  SHEATH POTENTIAL -E*DPHI (EV) FOR AN NP COMPONENT
C  PLASMA AS FUNCTION OF:
C  TE ELECTRON TEMPERATUR (EV)
C  DP(J),J=1,NP  BULK ION DENSITIES (#/CM**3)
C  VP(J),J=1,NP  BULK ION DRIFT VELOCITIES (CM/SEC)
C  NZP(J),J+1,NP CHARGE NUMBERS OF BULK IONS
C  GAMMA         SECONDARY ELECTRON EMISSION COEF.
C                SHEATH VOLTAGE
C                "SHEATH" DECREASES WITH INCREASING GAMMA
C  CUR           NET CURRENT DENSITY OF BULK PLASMA TO TARGET
C                CUR=(J-ION) - (J-ELECTR.) (AMP/CM**2)
C                SHEATH DECREASES FOR NEGATIVE CURRENTS, AND
C                       INCREASES FOR POSITIVE CURRENTS
C
      INCLUDE 'CCONA'
      DIMENSION DP(*),VP(*),NZP(*)
      SAVE ICOUNT
      DATA ICOUNT/0/
C
      SHEATH=0.
      DE=0.
      DO 10 J=1,NP
        DE=DE+DP(J)*NZP(J)
10    CONTINUE
      DE=DE+EPS60
C  UNITS OF SUM: VELOCITY (CM/SEC)
      SUM=0.
      DO 100 J=1,NP
        SUM=SUM+NZP(J)*DP(J)/DE*VP(J)
100   CONTINUE
      CE=CVEL2A*SQRT(TE/PMASSE)
      SUM=1./CE*SQRT(PI2A)/(1.-GAMMA)*(SUM-CUR/ELCHA/DE)
      IF (SUM.GT.0.D0) THEN
        SHEATH=-TE*LOG(SUM)
      ELSEIF (ICOUNT.LE.10) THEN
        WRITE (6,*) 'WARNING FROM FCT. SHEATH: INVALID ARGUMENTS '
        WRITE (6,*) 'ZERO SHEATH POTENTIAL RETURNED FOR SURFACE ', MS
        ICOUNT=ICOUNT+1
      ENDIF
      RETURN
      END
C
      SUBROUTINE VELOCS (TIWL,ESHET,VWL,VXWL,VYWL,VZWL,RSQDV,CVRSS,
     .                   CX,CY,CZ,
     .                   E0S,VELXS,VELYS,VELZS,VELS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  FETCH A NEW VELOCITY FROM A MAXWELLIAN FLUX AT A SURFACE GIVEN
C  BY THE NORMAL: CRTX,CRTY,CRTZ
C  THE MAXWELLIAN IS DEFINED BY A TEMPERATURE TIWL (EV) AND A
C  DRIFT VECTOR VXWL,VYWL,VZWL (CM/SEC). VWL IS THE VECTOR
C  NORM OF THIS DRIFT VECTOR, AND IS NEEDED ONLY TO DECIDE, WHETHER
C  THERE IS A DRIFT AT ALL (VWL.GT.0) OR NOT (VWL.LE.0.).
C  THE PARTICLE MAY BE ACCELERATED PERPENDICULAR
C  TOWARDS THE TARGET BY A SHEATH POTENTIAL WITH ENERGY ESHET (EV)
C  (IF ESHET.GT.0.)
C  OUTPUT ENERGY, SPEED UNIT VECTOR AND VELOCITY ARE, RESP.:
C         E0S,    VELXS,VELYS,VELZS AND VELS
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CRAND'
      INCLUDE 'CCONA'
      INCLUDE 'COMPRT'
C---------------------------------------------------------------------
C
      IF (INIV1.EQ.0) CALL FMAXWL
C
      ZARG2=SQRT(TIWL)*RSQDV
      ZARG=ZARG2*SQ2I
      IF (VWL.GT.0.D0) THEN
        CALL ROTATI(VXWL,VYWL,VZWL,VXDR,VYDR,VZDR,CX,CY,CZ)
      ELSE
        VXDR=0.
        VYDR=0.
        VZDR=0.
      ENDIF
C
C  MAXWELLIAN FLUX + DRIFT CONTRIBUTION
      VLLX=FM1(INIV1)*ZARG
      VLLY=FM2(INIV1)*ZARG+VYDR
      VLLZ=FM3(INIV1)*ZARG+VZDR
      INIV1=INIV1-1
C  WEIGHT CORRECTION DUE TO DRIFT-COMPONENT IN MAXWELLIAN FLUX:
C  NORMALIZED DRIFT COMPONENT    VMX=VXDR/ZARG2
C  NORMALIZED THERMAL COMPONENT  VLX=VLLX/ZARG2
      IF (VXDR.NE.0.D0) THEN
        VMX=VXDR/ZARG2
        CCM=0.6026*VMX
        IF (CCM.GE.1.) THEN
C         WRITE (6,*) 'WARNING FROM SUBR. VELOCS:'
C         WRITE (6,*) 'MACH NUMBER PERP. TO TARGET TOO LARGE FOR'
C         WRITE (6,*) 'RANDOM SAMPLING ALGORITHM, M-PERP= ',VMX
C         WRITE (6,*) 'ARTIFICIAL CUT OFF IS USED! '
C         CCM=0.95
C         VMX=CCM/0.6026
CDR USE REJECTION TECHNIQUE RATHER THAN WEIGHT CORRECTION TECHNIQUE
C         WEIGHT=1.
          VLLX=VXDR
          SHIFT=VMX
100       ARBV=RANF_EIRENE()*(SHIFT+4.)
          A1=SQRT(SHIFT**2+2.)
          A2=SHIFT+A1
          A3=SHIFT-A1
          A4=0.5+0.5*SHIFT*A3
          A5=A4-(ARBV-SHIFT)**2
          A6=EXP(A5)
          VFKT=2.*ARBV/A2*A6
          IF(RANF_EIRENE().GT.VFKT) GOTO 100
          VLLX=ARBV*ZARG2
          GOTO 1000
CDR
        ENDIF
        RCCM=1./(1.-CCM)
        VLLX=VLLX*SQRT(RCCM)
        VLX=VLLX/ZARG2
C
        VMXSQ=-VMX*VMX
        FNOM=EXP(VLX*(VMX+VMX-CCM*VLX))*RCCM
        FACTOR=FNOM/(1.+VMX*PISQ*(1.+ERF(VMX))/EXP(VMXSQ))
        WEIGHT=WEIGHT*FACTOR
      ENDIF
1000  CONTINUE
C
C  SHEATH CONTRIBUTION
      IF (ESHET.GT.0.D0) THEN
        VELSH=SQRT(ESHET)*RSQDV
        VLLX=SQRT(VELSH**2+VLLX**2)
      ENDIF
      CALL ROTATF(VELXS,VELYS,VELZS,VLLX,VLLY,VLLZ,CX,CY,CZ)
C
      VELSQ=VELXS*VELXS+VELYS*VELYS+VELZS*VELZS
      VELS=SQRT(VELSQ)
      VELXS=VELXS/VELS
      VELYS=VELYS/VELS
      VELZS=VELZS/VELS
      E0S=CVRSS*VELSQ
C
      RETURN
      END
C
      SUBROUTINE VELOCX(K,VXO,VYO,VZO,VLO,IOLD,NOLD,VELQ,NFLAG,
     .                  IRCX,DUMT,DUMV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  NFLAG= 1:       SAMPLING FROM MONOENERGETIC DISTRIBUTION
C                  OF ION SPEED IN 1D, X DIRECTION
C                  (I.E., DOUBLE DELTA FUNCTION)
C                  E=ESIGCX(IRCX,1)
C  NFLAG= 1:       SAMPLING FROM MONOENERGETIC DISTRIBUTION
C                  OF ION SPEED IN 2D, X-Y DIRECTION
C                  ISOTROPIC IN 2D
C                  E=ESIGCX(IRCX,1)
C  NFLAG= 1:       SAMPLING FROM MONOENERGETIC DISTRIBUTION
C                  OF ION SPEED FUNCTION IN 3D
C                  ISOTROPIC IN 3D
C                  E=ESIGCX(IRCX,1)
C  NFLAG= 2:       SAMPLING FROM SHIFTED MAXWELLIAN
C                  "FMAXW" AT TI AND V-DRIFT IN CELL K
C  NFLAG= 3:       SAMPLING FROM SHIFTED MAXWELLIAN + WEIGHT CORRECTION
C                  FACTOR = SIGMA*VREL*FMAXW/<SIGMA*VREL>
C                  OR ALTERNATIVELY: REJECTION
C
C  USED E.G. FOR VOLUME RECOMBINATION SOURCE (NFLAG=2)
C  OR TO  FETCH A NEW  VELOCITY FOR A NEUTRAL ATOM "IATM",
C  A NEUTRAL MOLECULE "IMOL" OR A TEST ION "IION"
C  AFTER CX-EVENT WITH BULK ION "IPLS" IN CELL NO. K FROM A SHIFTED
C  MAXWELLIAN (NFLAG=2), WEIGHTED BY SIGMA*VREL (NFLAG=3)
C  K   : .NE.0 :CELL INDEX FOR LOCAL BULK ION TI AND V_DRIFT
C        .EQ.0 :TX,TY,TZ,V-DRIFT_X,Y,Z ARE NOT FROM LOCAL BULK ION
C               SPECIES IPLS, BUT EXPLICITLY IN THE PARAMETERS DUMT AND DUMV.
C               RESPECTIVELY.
C  VXO : X COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VYO : Y COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VZO : Z COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VLO : VELOCITY OF TEST PARTICLE BEFORE EVENT
C  IOLD: SPECIES INDEX OF THE TEST PARTICLE BEFORE THE EVENT
C  NOLD: DITO, IN MODCOL-ARRAY
C  IPLS: SPECIES INDEX FOR THE THERMAL PLASMA ION VELOCITY
C        AND FOR THE PLASMA DRIFT VELOCITY TO BE USED AS
C        SHIFT VECTOR   (IPLS IN COMMON COMUSR)
C  IRCX: LABEL FOR REACTION, E.G., FOR SIGVCX(IRCX)
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CRAND'
      INCLUDE 'COMXS'
      INCLUDE 'CZT1'
      INCLUDE 'COMUSR'
      INCLUDE 'CLOGAU'
      INCLUDE 'CCONA'
      DIMENSION DUMT(3),DUMV(3)
C
      IF (K.GT.0) THEN
        ZARGX=ZRG(IPLS,K)
        ZARGY=ZRG(IPLS,K)
        ZARGZ=ZRG(IPLS,K)
        IF (NLDRFT) THEN
          VXDR=VXIN(IPLS,K)
          VYDR=VYIN(IPLS,K)
          VZDR=VZIN(IPLS,K)
        ELSE
          VXDR=0.D0
          VYDR=0.D0
          VZDR=0.D0
        ENDIF
      ELSE
        IF (NFLAG.NE.2) GOTO 999
        ZARGX=DUMT(1)
        ZARGY=DUMT(2)
        ZARGZ=DUMT(3)
        VXDR=DUMV(1)
        VYDR=DUMV(2)
        VZDR=DUMV(3)
      ENDIF
C
      IF (INIV2.LE.0) CALL FGAUSS
C
C  SET NEW VELOCITY
      VXN=FG1(INIV2)*ZARGX+VXDR
      VYN=FG2(INIV2)*ZARGY+VYDR
      VZN=FG3(INIV2)*ZARGZ+VZDR
      INIV2=INIV2-1
C
C
C   DEFAULT MODEL
C
      IF (NFLAG.EQ.2) THEN
C
        VELQ=VXN*VXN+VYN*VYN+VZN*VZN
        VEL=SQRT(VELQ)
        VN=1./VEL
        VELX=VXN*VN
        VELY=VYN*VN
        VELZ=VZN*VN
C
C       VEL_MEAN=SQRT(VXDR**2+VYDR**2+VZDR**2)
C       VELX_MEAN=VXDR/(VEL_MEAN+EPS60)
C       VELY_MEAN=VYDR/(VEL_MEAN+EPS60)
C       VELZ_MEAN=VZDR/(VEL_MEAN+EPS60)
C
C   NOTHING MORE TO BE DONE
C
      ELSEIF (NFLAG.EQ.3) THEN
C
C   WEIGHT CORRECTION DUE TO ENERGY DEPENDENCE IN CROSS SECTION
C   (IF DATA ARE AVAILABLE)
C
        VX=VXO*VLO
        VY=VYO*VLO
        VZ=VZO*VLO
        VRELQ=(VXN-VX)**2+(VYN-VY)**2+(VZN-VZ)**2
        VREL=SQRT(VRELQ)
C
        ELAB=LOG(VRELQ)+DEFCX(IRCX)
        IREAC=MODCOL(3,1,NOLD,IPLS)
        CXS=CROSS(ELAB,IREAC)
        WEIGHT=WEIGHT*CXS*VREL*DIIN(IPLS,K)/SIGVCX(IRCX)
C
        VELQ=VXN*VXN+VYN*VYN+VZN*VZN
        VEL=SQRT(VELQ)
        VN=1./VEL
        VELX=VXN*VN
        VELY=VYN*VN
        VELZ=VZN*VN
C
C       VEL_MEAN=VEL
C       VELX_MEAN=VELX
C       VELY_MEAN=VELY
C       VELZ_MEAN=VELZ
C
C   1D THERMAL ION SPEED (FOR 1D 1 ENERGY GROUP APPROXIMATION)
C
      ELSEIF (NFLAG.EQ.1) THEN
        write (6,*) 'warning, nflag=1 in velocx '
        VEL=RSQDVP(IPLS)*SQRT(ESIGCX(IRCX,1))
        VELQ=VEL*VEL
        VN=1./SQRT(VXN*VXN)
        VELX=VXN*VN
        VELY=0.
        VELZ=0.
C
C       VEL_MEAN=VEL
C       VELX_MEAN=VELX
C       VELY_MEAN=VELY
C       VELZ_MEAN=VELZ
C
C   2D THERMAL ION SPEED (FOR 2D 1 ENERGY GROUP APPROXIMATION)
C
C     ELSEIF (NFLAG.EQ.1) THEN
C       VEL=RSQDVP(IPLS)*SQRT(ESIGCX(IRCX,1))
C       VELQ=VEL*VEL
C       VN=1./SQRT(VXN*VXN+VYN*VYN)
C       VELX=VXN*VN
C       VELY=VYN*VN
C       VELZ=0.
C
C   RESET TO 3D THERMAL ION SPEED (FOR 3D 1 ENERGY GROUP APPROXIMATION)
C
C     ELSEIF (NFLAG.EQ.1) THEN
C       VEL=RSQDVP(IPLS)*SQRT(ESIGCX(IRCX,1))
C       VELQ=VEL*VEL
C       VN=1./SQRT(VXN*VXN+VYN*VYN+VZN*VZN)
C       VELX=VXN*VN
C       VELY=VYN*VN
C       VELZ=VZN*VN
      ENDIF
C
      RETURN
C
999   CONTINUE
      WRITE (6,*) 'PARAMETER ERROR IN SUBR. VELOCX. EXIT CALLED'
      CALL EXIT
      END
C
      SUBROUTINE VELOEI(K,IRDS,VXO,VYO,VZO,VLO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  FETCH A NEW VELOCITY OF TEST PARTICLE AFTER ELECTRON IMPACT COLLISION
C
C  K   : CELL INDEX
C  VXO : X COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VYO : Y COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VZO : Z COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VLO : VELOCITY OF TEST PARTICLE BEFORE EVENT
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMUSR'
      INCLUDE 'COMXS'
      INCLUDE 'CRAND'
      INCLUDE 'CZT1'
      INCLUDE 'CCONA'
C
C  FIND TYPE OF NEXT GENERATION PARTICLE (ATOM, MOLECULE, TEST ION)
C
      ZEP3=RANF_EIRENE( )
C
c slmod begin - not tr
c      IF (output)
c     .WRITE(6,'(A,2I6,3(I6,F12.6))')
c     .  'MARK: P2ND ',K,IRDS,
c     .  NSPA  ,P2ND(IRDS,NSPA  ),
c     .  NSPAM ,P2ND(IRDS,NSPAM ),
c     .  NSPAMI,P2ND(IRDS,NSPAMI)
c slmod end
      IF (ZEP3.LE.P2ND(IRDS,NSPA)) THEN
C
C  A NEUTRAL ATOM IS BORN, FIND SPECIES INDEX IATM AND WEIGHT
C
        ITYP=1
        DO 448 IATM=1,NATMIM
          ISPZA=IATM
          IF (ZEP3.LE.P2ND(IRDS,ISPZA)) GOTO 449
448     CONTINUE
        IATM=NATMI
449     CONTINUE
        CVRSS=CVRSSA(IATM)
        RSQDV=RSQDVA(IATM)
        EDISS=EATDS(IRDS,IATM,2)
C
      ELSEIF (ZEP3.LE.P2ND(IRDS,NSPAM)) THEN
C
C  A NEUTRAL MOLECULE IS BORN, FIND SPECIES INDEX IMOL AND WEIGHT
C
        ITYP=2
        DO 458 IMOL=1,NMOLIM
          ISPZM=NATMI+IMOL
          IF (ZEP3.LE.P2ND(IRDS,ISPZM)) GOTO 459
458     CONTINUE
        IMOL=NMOLI
459     CONTINUE
        CVRSS=CVRSSM(IMOL)
        RSQDV=RSQDVM(IMOL)
        EDISS=EMLDS(IRDS,IMOL,2)
C
      ELSEIF (ZEP3.LE.P2ND(IRDS,NSPAMI)) THEN
C
C  A TEST ION IS BORN, FIND SPECIES INDEX IION AND WEIGHT
C
        ITYP=3
        DO 468 IION=1,NIONIM
          ISPZI=NSPAM+IION
          IF (ZEP3.LE.P2ND(IRDS,ISPZI)) GOTO 469
468     CONTINUE
        IION=NIONI
469     CONTINUE
        CVRSS=CVRSSI(IION)
        RSQDV=RSQDVI(IION)
        EDISS=EIODS(IRDS,IION,2)
C
      ELSE
        WRITE (6,*) 'ERROR IN VELOEI '
        WRITE (6,*) 'IREI ',IRDS,P2ND(IRDS,NSPAMI)
        CALL EXIT
      ENDIF
C
      IF (NSTORDR >= NRAD) THEN
        EHEAVY=EHVDS1(IRDS,K)
      ELSE
        EHEAVY=FEHVDS1(IRDS,K)
      END IF
      EDISS=EDISS*EHEAVY
C
C  FIND SPEED VECTOR FROM ISOTROPIC DISTRIBUTION IN CENTER OF MASS
C  SYSTEM
C

      IF (EDISS.GT.0.D0) THEN
C
C  NEXT LINES: E-NEW=E-OLD+EDIS, ON AVERAGE
C
        IF (INIV3.EQ.0) CALL FISOTR
C
        VXDIS=FI1(INIV3)
        VYDIS=FI2(INIV3)
        VZDIS=FI3(INIV3)
        INIV3=INIV3-1
C
        VELDS=SQRT(EDISS)*RSQDV
        VX=VLO*VXO+VELDS*VXDIS
        VY=VLO*VYO+VELDS*VYDIS
        VZ=VLO*VZO+VELDS*VZDIS
        VELQ=VX*VX+VY*VY+VZ*VZ
        VEL=SQRT(VELQ)
        VELX=VX/VEL
        VELY=VY/VEL
        VELZ=VZ/VEL
      ELSE
        VELX=VXO
        VELY=VYO
        VELZ=VZO
        VEL=VLO
        VELQ=VEL*VEL
      ENDIF
      E0=CVRSS*VELQ
C
      RETURN
      END
C
C
      SUBROUTINE VELOEL(K,VXO,VYO,VZO,VLO,IOLD,NOLD,VELQ,NFLAG,
     .                  IREL,RMASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  THIS SUBROUTINE CARRIES OUT AN ELASTIC COLLISION OF A TEST PARTICLE
C  WITH A BULK PARTICLE.
C  IT RETURNS THE POST COLLISION VELOCITY VECTOR.
C
C  1ST STEP: FIND COLLISION PARTNER FROM BULK ION SPECIES "IPLS":
C            (VXN,VYN,VZN)
C  NFLAG= 1:       SAMPLING FROM MONOENERGETIC DISTRIBUTION
C                  OF ION SPEED IN 1D, X DIRECTION
C                  (I.E., DOUBLE DELTA FUNCTION)
C                  E=ESIGCX(IRCX,1)
C  NFLAG= 1:       SAMPLING FROM MONOENERGETIC DISTRIBUTION
C                  OF ION SPEED IN 2D, X-Y DIRECTION
C                  ISOTROPIC IN 2D
C                  E=ESIGCX(IRCX,1)
C  NFLAG= 1:       SAMPLING FROM MONOENERGETIC DISTRIBUTION
C                  OF ION SPEED FUNCTION IN 3D
C                  ISOTROPIC IN 3D
C                  E=ESIGCX(IRCX,1)
C  NFLAG= 2:       SAMPLING FROM SHIFTED MAXWELLIAN
C                  "FMAXW" AT TI AND V-DRIFT IN CELL K
C  NFLAG= 3:       SAMPLING FROM SHIFTED MAXWELLIAN + WEIGHT CORRECTION
C                  FACTOR = SIGMA*VREL*FMAXW/<SIGMA*VREL>
C                  OR ALTERNATIVELY: REJECTION
C
C  2RD STEP: FIND IMPACT PARAMETER B AND/OR SCATTERING ANGLE
C  3TH STEP: FIND NEW VELOCITY VECTOR
C
C
C  K   : CELL INDEX
C  VXO : X COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VYO : Y COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VZO : Z COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VLO : VELOCITY OF TEST PARTICLE BEFORE EVENT
C  IOLD: SPECIES INDEX OF THE TEST PARTICLE BEFORE THE EVENT
C  NOLD: DITO, IN MODCOL-ARRAY
C  IPLS: SPECIES INDEX FOR PLASMA DRIFT VELOCITY TO BE USED AS
C        SHIFT VECTOR AND FOR THERMAL ION SPEED DEFINITION
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CRAND'
      INCLUDE 'COMXS'
      INCLUDE 'CZT1'
      INCLUDE 'COMUSR'
      INCLUDE 'CLOGAU'
      INCLUDE 'CCONA'
      INCLUDE 'CTRCEI'
      COMMON /PAR/ B,ER,IFLAG
      LOGICAL LGWGHT,LLAST
      DIMENSION XCMEAN(NREL),NCMEAN(NREL),IFLREL(NREL)
      data wrmax/-1000.d0/,wrmin/1000.d0/ipmax/0/,ipmin/0/
      DATA IFIRST/0/
      SAVE
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        DO IRL=1,NRELI
          IFLREL(IRL)=0
          XCMEAN(IRL)=0.D0
          NCMEAN(IRL)=0
        ENDDO
        LLAST=.FALSE.
      ENDIF
C
      LLAST=LLAST.AND.LGLAST
C
      IF (IFLREL(IREL).EQ.0) THEN
        IFLREL(IREL)=-1
C  PREPARE REJECTION SAMPLING OF INCIDENT ION VELOCITY
C  IS CROSS SECTION AVAILABLE?
        IREAC=MODCOL(5,1,NOLD,IPLS)
        IF (IREAC.EQ.0) GOTO 1
C
        elmin=log(0.01)
        elmax=log(1.e3)
        SIGVMX=-1.D60
        JJ=1
        do j=1,1000
          elab=elmin+(j-1)/999.*(elmax-elmin)
          CEL=CROSS(ELAB,IREAC)
          vrq=exp(elab-defel(IREL))
          vr=sqrt(vrq)
          if (cel*vr.gt.SIGVMX) then
            JJ=J
            SIGVMX=cel*vr
          endif
        enddo
        WRITE (6,*) 'IREL,SIGVMX IN VELOEL,JJ ',IREL,SIGVMX,JJ
        IF (JJ.NE.1.AND.JJ.NE.1000) IFLREL(IREL)=1
      ENDIF
C
1     CONTINUE
c
      ICOUNT=1
123   CONTINUE
      ZARG=ZRG(IPLS,K)
      VX=VXO*VLO
      VY=VYO*VLO
      VZ=VZO*VLO
C
      IF (INIV2.LE.0) CALL FGAUSS
C
C  NEXT: STEP 1
C
      IF (NLDRFT) THEN
        VXDR=VXIN(IPLS,K)
        VYDR=VYIN(IPLS,K)
        VZDR=VZIN(IPLS,K)
      ELSE
        VXDR=0.D0
        VYDR=0.D0
        VZDR=0.D0
      ENDIF
C
      VXN=FG1(INIV2)*ZARG+VXDR
      VYN=FG2(INIV2)*ZARG+VYDR
      VZN=FG3(INIV2)*ZARG+VZDR
      INIV2=INIV2-1
C
C  SAMPLE FROM MAXWELLIAN OR FROM WEIGHTED MAXWELLIAN?
      LGWGHT=IFLREL(IREL).GT.0.AND.NFLAG.EQ.3
C
      IF (.NOT.LGWGHT) THEN
C
        VXI=VXN
        VYI=VYN
        VZI=VZN
C
      ELSEIF (LGWGHT) THEN
C
C   WEIGHT CORRECTION DUE TO ENERGY DEPENDENCE IN CROSS SECTION
C   OR: REJECTION     DUE TO ENERGY DEPENDENCE IN CROSS SECTION
C   PRESENT VERSION: REJECTION
        VRELQ=(VXN-VX)**2+(VYN-VY)**2+(VZN-VZ)**2
        VREL=SQRT(VRELQ)
        ELAB=LOG(VRELQ)+DEFEL(IREL)
        IREAC=MODCOL(5,1,NOLD,IPLS)
        CEL=CROSS(ELAB,IREAC)
C
c       if (nlrejc) then
          TEST=RANF_EIRENE()*SIGVMX
          if (test.gt.cel*vrel) then
c  reject
            icount=icount+1
            if (icount.lt.500) goto 123
            write (6,*) 'icount too large IN VELOEL. ACCEPT SAMPLE '
          else
c  accept
            xcmean(irel)=xcmean(irel)+icount
            ncmean(irel)=ncmean(irel)+1
            if (lglast.and..not.llast.and.TRCLST) THEN
              llast=.true.
              write (6,*) 'sampling efficiancy in veloel '
              write (6,*) 'IREL, MEAN NUMBER OF SAMPLINGS'
              do irel=1,nreli
                write (6,*) irel,xcmean(irel)/(ncmean(irel)+eps60)
              enddo
            endif
          endif
c       elseif (nlweight) then
C         wo=weight
c         WEIGHT=WEIGHT*CEL*VREL*DIIN(IPLS,K)/SIGVEL(IREL)
c
c         wrat=wo/weight
C         wrmean=wrmean+wrat
c         imean=imean+1
c         if (wrat.gt.wrmax) then
c           ipmax=npanu
c           wrmax=wrat
c           timax=tiin(ipls,k)
c           e0max=e0
c           elabs=elab
c           vrels=vrel
c           sigs=sigvel(IREL)/diin(ipls,k)
c         endif
c         if (wrat.lt.wrmin) then
c           ipmin=npanu
c           wrmin=wrat
c         endif
c         if (lglast) write (6,*) 'wrmax... ',wrmax,ipmax,wrmin,ipmin
c         if (lglast) write (6,*) 'timax... ',timax,e0max
c         if (lglast) write (6,*) 'elabs... ',elabs,vrels,sigs
c         if (lglast) write (6,*) 'mean... ',wrmean/(float(imean)+1.d-50),
c    .                                       xcmean/(float(ncmean)+1.d-50)
c       endif
C
        VXI=VXN
        VYI=VYN
        VZI=VZN
C
C   1D THERMAL ION SPEED (FOR 1D 1 ENERGY GROUP APPROXIMATION)
C
      ELSEIF (NFLAG.EQ.1) THEN
        write (6,*) 'warning: nflag=1 in veloel'
        VL=ZARG
        VN=VL/SQRT(VXN*VXN)
        VXI=VXN*VN
        VYI=0.
        VZI=0.
C
C   2D THERMAL ION SPEED (FOR 2D 1 ENERGY GROUP APPROXIMATION)
C
C     ELSEIF (NFLAG.EQ.1) THEN
C       VL=SQRT(2.)*ZARG
C       VN=VL/SQRT(VXN*VXN+VYN*VYN)
C       VXI=VXN*VN
C       VYI=VYN*VN
C       VZI=0.
C
C   3D THERMAL ION SPEED (FOR 3D 1 ENERGY GROUP APPROXIMATION)
C
C     ELSEIF (NFLAG.EQ.1) THEN
C       VL=SQRT(3.)*ZARG
C       VN=VL/SQRT(VXN*VXN+VYN*VYN+VZN*VZN)
C       VXI=VXN*VN
C       VYI=VYN*VN
C       VZI=VZN*VN
      ENDIF
C
C  STEP 1 FINISHED, COLLIDING BULK ION'S VELOCITY IS SET: VXI,VYI,VZI
C
200   CONTINUE
C
C  FIND TYPE OF COLLISION: IFLAG
C
C  NEUTRAL ATOMS:
C
      IF (ITYP.EQ.1) THEN
C  NEUTRAL-NEUTRAL, IN BGK APPROXIMATION
        IF (NCHRGP(IPLS).EQ.0) THEN
          IFLAG=0
C  H+ + H, IFLAG=1 ?
        ELSEIF (NCHARA(IATM).EQ.1.AND.NCHARP(IPLS).EQ.1.AND.
     .      NCHRGP(IPLS).EQ.1) THEN
          IFLAG=1
C  H+ + HE, IFLAG=2 ?
        ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.1.AND.
     .      NCHRGP(IPLS).EQ.1) THEN
          IFLAG=2
C  HE+ + HE, IFLAG=4 ?
        ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.2.AND.
     .      NCHRGP(IPLS).EQ.1) THEN
          IFLAG=4
        ELSE
          GOTO 995
        ENDIF
C
C  NEUTRAL MOLECULES:
C
      ELSEIF (ITYP.EQ.2) THEN
C  NEUTRAL-NEUTRAL, BGK APPROXIMATION
        IF (NCHRGP(IPLS).EQ.0) THEN
          IFLAG=0
C  H+ + H2, IFLAG=3 ?
        ELSEIF (NCHARM(IMOL).EQ.2.AND.NCHARP(IPLS).EQ.1.AND.
     .      NCHRGP(IPLS).EQ.1) THEN
          IFLAG=3
        ELSE
          GOTO 995
        ENDIF
      ELSE
        GOTO 995
      ENDIF
C
C  NEXT: STEP 2, FIND PRE COLLISION DATA AND TOTAL CROSS SECTION
C
      IF (IFLAG.NE.0) THEN
        RMN=RMASS
        RMI=RMASSP(IPLS)
        RMSI=1./(RMN+RMI)
        RLMS=RMN*RMI*RMSI
C  RELATIV VELOCITY AND RELATED DATA
        VRELX=VX-VXI
        VRELY=VY-VYI
        VRELZ=VZ-VZI
        VRQYZ=         VRELY**2+VRELZ**2
        VRELQ=VRELX**2+VRQYZ
        ER=RLMS*VRELQ*CVELI2
        VREL=SQRT(VRELQ)
        VRYZ=SQRT(VRQYZ+EPS60)
C  CENTER OF MASS VELOCITY
        VSX=(RMI*VXI+RMN*VX)*RMSI
        VSY=(RMI*VYI+RMN*VY)*RMSI
        VSZ=(RMI*VZI+RMN*VZ)*RMSI
C  TOTAL CROSS SECTION (ONLY IF NOT COMPUTED EARLIER AT THIS CALL)
        IF (.NOT.LGWGHT) THEN
          ELAB=LOG(VRELQ)+DEFEL(IREL)
          IREAC=MODCOL(5,1,NOLD,IPLS)
          IF (IREAC.EQ.0) GOTO 995
          CEL=CROSS(ELAB,IREAC)
        ENDIF
C
C  STEP 2 FINISHED, CROSS SECTION CEL IS FOUND
C  NEXT: STEP 3
C
        BMAX=SQRT(CEL*PIAI)/0.52917E-8
        B= SQRT(RANF_EIRENE( ))*BMAX
      ENDIF
C
C  STEP 3 FINISHED, IMPACT PARAMETER IS FOUND , IN UNITS: BOHR RADIA
C  NEXT: STEP 4
C
      IF (IFLAG.EQ.0) THEN
C
C  THIS PART: ONLY RELAXATION TO MAXWELLIAN, I.E., POST COLLISION
C             NEUTRAL SAMPLED FROM BULK ION POPULATION
C
        VELQ=VXI*VXI+VYI*VYI+VZI*VZI
        VEL=SQRT(VELQ)
        VELX=VXI/VEL
        VELY=VYI/VEL
        VELZ=VZI/VEL
C
C       IF (.NOT.LGWGHT) THEN
C         VEL_MEAN=SQRT(VXDR**2+VYDR**2+VZDR**2)
C         VELX_MEAN=VXDR/(VEL_MEAN+EPS60)
C         VELY_MEAN=VYDR/(VEL_MEAN+EPS60)
C         VELZ_MEAN=VZDR/(VEL_MEAN+EPS60)
C       ELSE
C         VEL_MEAN=VEL
C         VELX_MEAN=VELX
C         VELY_MEAN=VELY
C         VELZ_MEAN=VELZ
C       ENDIF

      ELSE
C
C  COLLISION PARAMETERS IFLAG,ER AND B ARE DEFINED NOW.
C
C  FIND DISTANCE OF CLOSEST APPROACH: RSTERN
C
        RS=RSTERN(ER,B,IFLAG)
C
C  INTEGRAL TO FIND DEFLECTION ANGLE CHI
C
        CALL GAUMEH (RS,ER,B,IFLAG,10,1,RESULT)
        CHI=PIA-2.*B/RS*RESULT
C
C  CONVERT FROM DEFLECTION ANGLE TO OBSERVABLE SCATTERING ANGLE
        PH=ACOS(COS(CHI))
C
C  POLAR ANGLE
        EPS=PI2A*RANF_EIRENE( )
C
        CPH=COS(PH)
        SPH=SIN(PH)
        CEPS=COS(EPS)
        SEPS=SIN(EPS)
        VRSX=VRELX*CPH+SPH*SEPS*VRYZ
        VRSY=VRELY*CPH+SPH*(VREL*VRELZ*CEPS-VRELX*VRELY*SEPS)/VRYZ
        VRSZ=VRELZ*CPH-SPH*(VREL*VRELY*CEPS+VRELX*VRELZ*SEPS)/VRYZ
C
        VELX=VSX+RLMS/RMN*VRSX
        VELY=VSY+RLMS/RMN*VRSY
        VELZ=VSZ+RLMS/RMN*VRSZ
        VELQ=VELX*VELX+VELY*VELY+VELZ*VELZ
        VEL=SQRT(VELQ)
        VELX=VELX/VEL
        VELY=VELY/VEL
        VELZ=VELZ/VEL
C
C       VEL_MEAN=VEL
C       VELX_MEAN=VELX
C       VELY_MEAN=VELY
C       VELZ_MEAN=VELZ
      ENDIF
C
C  STEP 4 FINISHED, POST COLLISION VELOCITY IS SET
C  NEXT: RETURN
      IF (output) WRITE(0,*) 'MARK: RETURN FROM VELOEL'
      RETURN
C
995   CONTINUE
      WRITE (6,*) 'ERROR IN VELOEL, NO ELASTIC COLLISION DATA'
      WRITE (6,*) 'AVAILABLE '
      WRITE (6,*) 'ITYP,IATM,IMOL,IION,IPLS ',ITYP,IATM,IMOL,IION,IPLS
      CALL EXIT
      END
C
      SUBROUTINE FOLNEUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     NEUTRAL PARTICLE, LAUNCHED AT X0,Y0,Z0, IN CELL NRCELL, IPOLG,
C     NPCELL, NTCELL, NACELL, NBLOCK, IS FOLLOWED
C
C  ON INPUT
C     ITYP=1 OR ITYP=2
C  ON OUTPUT:
C
C     LGPART=TRUE
C     ITYP=3        IF A NEXT GENERATION TEST ION IS BORN
C                   IION= SPECIES INDEX OF NEXT GENERATION ION
C     ITYP=0 OR ITYP=4 AND LGPART=FALSE  ELSE
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'COMNNL'
      INCLUDE 'CZT1'
      INCLUDE 'CESTIM'
      INCLUDE 'CGRID'
c slmod begin - not tr (checking new version to see if required)
      INCLUDE 'CGEOM'
      INCLUDE 'CPOLYG'
c slmod end
      INCLUDE 'CCONA'
      INCLUDE 'COMUSR'
      INCLUDE 'COUTAU'
      INCLUDE 'CLOGAU'
      INCLUDE 'CTRIG'
      INCLUDE 'CLGIN'
      INCLUDE 'COMXS'
      INCLUDE 'CSPEZ'
      INCLUDE 'COMSPL'
      INCLUDE 'CRAND'
c slmod begin - debug - tr
      INTEGER PolPos,CheckCell,ip1,ip2,ii
      DIMENSION cp(4)

      COMMON /TRASHCOM/ nlost
      INTEGER           nlost

      COMMON /STATSCOM/ npar,xpos,ypos,zpos,xvel,yvel,zvel
      INTEGER           npar(9)
      REAL              xpos(9),ypos(9),zpos(9),
     .                  xvel(9),yvel(9),zvel(9)

      INTEGER npanu_last,slcount
      REAL*8 partime(2)

      SAVE

      DATA partime /0.0D0,0.0D0/
c slmod end
      DIMENSION CFLAG(6,3)
      DIMENSION AX(2),XSTOR(NSTOR),RPST(NPARTC),IPST(MPARTC)
      DIMENSION XSTOR2(NSTOR,N2ND+N3RD),XSTORC(NSTOR)
      LOGICAL NLPR,LTRANS
      EQUIVALENCE (XSTOR(1),SIGVCX(1))
C
C  TENTATIVELY ASSUME: A NEXT GENERATION PARTICLE WILL BE BORN
      LGPART=.TRUE.
      XGENER=0
C
C  THE  CELL NUMBER NRCELL, IPOLG, NPCELL, NTCELL, NACELL, NBLOCK
C  WAS ALREADY SET IN CALLING SUBROUTINE MCARLO
C
C  IF NLSRFX, SURFACE INDEX MRSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFY, SURFACE INDEX MPSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFZ, SURFACE INDEX MTSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFA, SURFACE INDEX MASURF MUST BE DEFINED AT THIS POINT
C
C  EACH NEUTRAL PARTICLE TRACK STARTS AT THIS POINT
C
c slmod begin - tr
 
c      WRITE(0,*) 'npanu',npanu,npanu_last

      IF (npanu.NE.npanu_last) THEN
        npanu_last = npanu
        npar(istra) = npar(istra) + 1
        xpos(istra) = xpos(istra) + x0
        ypos(istra) = ypos(istra) + y0
        zpos(istra) = zpos(istra) + z0
        xvel(istra) = xvel(istra) + velx
        yvel(istra) = yvel(istra) + vely
        zvel(istra) = zvel(istra) + velz
       
        slcount = 0


        CALL CLOCK(partime(1))

c        WRITE(0,*) 'setting partime1',partime(1)

      ENDIF

c      write(0,*) 'z0=',z0

c      WRITE(6,'(I5,1X,3F10.5,1X,1P,3E12.4,0P,1X,F10.5)')
c     .  npanu,x0,y0,z0,velx,vely,velz,weight
c slmod end
101   CONTINUE
c slmod begin - debug - tr

      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A,F12.4)') 'FOLNEUT: C E0=',e0



      IF (debugopt.GT.20) THEN
        IF (npanu.GE.debugopt.AND.npanu.LT.debugopt+100) THEN
          printopt = 1
        ELSE
          printopt = 0
        ENDIF
      ENDIF

c      IF (istra.EQ.3) THEN
c        printopt = 1
c      ENDIF

      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(3X,A)') 'FOLNEUT: '
        WRITE(6,'(3X,A,I7,F7.4)')
     .    'FOLNEUT: START OF TRACK (NPANU WEIGHT) ',
     .    npanu,weight
      ENDIF

      IF (debugopt.EQ.-3)
     .  WRITE(6,'(A,I6)') 'TRACK = ',npanu
c slmod end
      IF (ITYP.EQ.1) THEN
        IF (IATM.LE.0.OR.IATM.GT.NATMI) GOTO 998
        LOGATM(IATM,ISTRA)=.TRUE.
        NLPR=NLPRCA(IATM)
        NRC=NRCA(IATM)
      ELSEIF (ITYP.EQ.2) THEN
        IF (IMOL.LE.0.OR.IMOL.GT.NMOLI) GOTO 998
        LOGMOL(IMOL,ISTRA)=.TRUE.
        NLPR=NLPRCM(IMOL)
        NRC=NRCM(IMOL)
      ENDIF
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
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A,F12.4)') 'FOLNEUT: D E0=',e0
c slmod end
      NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
      NJUMP=0
      DO 102 J=1,NR1ST
102     TIMINT(J)=0.0
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
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(3X,A)') 'FOLNEUT: '
        WRITE(6,'(3X,2A)')
     .    'FOLNEUT: NRCELL,NPCELL,NTCELL,NP2T3,NR1P2,',
     .    'NBLCKA,NCELL,NLIMII,NLIMIE = '
        WRITE(6,'(7X,9I5)')
     .    nrcell,npcell,ntcell,np2t3,nr1p2,nblcka,ncell,
     .    nlimii(ncell),nlimie(ncell)
      ENDIF
c slmod end
C
C TL: DISTANCE TO NEXT ADDITIONAL SURFACE
c slmod begin - debug (temp) - not tr
c      IF (NACELL.GT.0) THEN
c        printopt = 1
c      ELSE
c        printopt = 0
c      ENDIF

c slmod end
      IF (NLIMII(NCELL).LE.NLIMIE(NCELL)) THEN
c slmod begin - debug - tr
        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(3X,A)') 'FOLNEUT: '
          WRITE(6,'(3X,A)') 'FOLNEUT: Calling TIMEA1'
          WRITE(6,'(3X,2A)')
     .      'FOLNEUT: Output (NR,NP IR,IP MRSURF ZDT1 ',
     .      'X0,Y0,Z0 X00,Y00,Z00)'
          WRITE(6,'(5X,A,2I4,A,2I4,A,I4,E11.3,2X,2(3F9.2,2X))')
     .      '(',NRCELL,NPCELL,') (',IRCELL,IPCELL,')',MRSURF,ZDT1,
     .      X0,Y0,Z0,X00,Y00,Z00
        ENDIF
c slmod end
        CALL TIMEA1 (MSURF,NCELL,NTCELL,X0,Y0,Z0,TIME,
     .               VELX,VELY,VELZ,VEL,
c slmod begin - tr - new
     .               MASURF,XLI,YLI,ZLI,SG,TL,NLTRC,NACELL,NTRSEG)
c
c     .               MASURF,XLI,YLI,ZLI,SG,TL,NLTRC)
c slmod end
        NLPR=NLPRCS(MASURF).OR.NLPR
        ZDT1=TL
        ZTST=TL
        CLPD(1)=ZDT1
        IF (MASURF.NE.0) ISRFCL=1
      ENDIF
C
C TT: DISTANCE UNTIL NEXT TIMESTEP LIMIT IS REACHED
      IF (LGTIME) THEN
        TT=(DTIMVI-TIME)*VEL
        IF (TT.LT.TL) THEN
          ZDT1=TT
          ZTST=TT
          CLPD(1)=ZDT1
          ISRFCL=2
        ENDIF
      ENDIF
C
C  SCAN OVER SEGMENT
C
210   CONTINUE
C
C  TS:   DISTANCE TO NEXT SURFACE OF STANDARD MESH
C  ZDT1: DISTANCE TRAVELLED IN CURRENT RADIAL CELL
C
c slmod begin - debug - tr

      IF (.FALSE.) THEN
        IF (nrcell.LT.1.OR.nrcell.GE.nr1st) THEN
          WRITE(6,'(3X,4I5,A)') NPANU,NPCELL,NRCELL,NACELL,' OUTSIDE'
c          WRITE(0,'(3X,4I5,A)') NPANU,NPCELL,NRCELL,NACELL,' OUTSIDE'
        ELSE
          WRITE(6,'(3X,4I5  )') NPANU,NPCELL,NRCELL,NACELL
c          WRITE(0,'(3X,4I5  )') NPANU,NPCELL,NRCELL,NACELL
        ENDIF
      ENDIF

      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(3X,A)') 'FOLNEUT: '
        WRITE(6,'(3X,A,15X,A)') 'FOLNEUT: ','*** START OF LOOP ***'
        WRITE(6,'(3X,A)') 'FOLNEUT: '
        IF (nrcell.LT.1.OR.nrcell.GE.nr1st) THEN
          WRITE(6,'(3X,A)') 'FOLNEUT: OUTSIDE GRID'
          WRITE(6,'(3X,A)') 'FOLNEUT: '
        ENDIF
        WRITE(6,'(3X,2A)')
     .    'FOLNEUT: Output (NR,NP NA IR,IP MRSURF ZDT1 ',
     .    'X0,Y0,Z0 X00,Y00,Z00)'
        WRITE(6,'(5X,A,2I4,A,I3,A,2I4,A,I4,E11.3,2X,2(3F9.2,2X))')
     .    '(',NRCELL,NPCELL,')',NACELL,' (',IRCELL,IPCELL,')',MRSURF,
     .    ZDT1,X0,Y0,Z0,X00,Y00,Z00
      ENDIF

      slcount = slcount + 1

      IF (slcount.EQ.1000) THEN
        CALL CLOCK(partime(2))     
c        WRITE(0,'(A,I7,1P,E10.2,0P)') 
c     .    'PARTIME:',npanu_last,partime(2)-partime(1)
        IF (partime(2)-partime(1).LT.600.0D0) THEN
          slcount = 0
        ELSE
c          WRITE(0,*) 'partime:',partime(1),partime(2)
          WRITE(6,*) 'TURNING ON DETAILED PARTICLE TRACKING'
          WRITE(0,*) 'TURNING ON DETAILED PARTICLE TRACKING'
          NLTRC=.TRUE.
        ENDIF
      ELSEIF (slcount.EQ.2000) THEN
        WRITE(0,*) 'EXITING EIRENE TO CHECK TRAPPED PARTICLE'
        WRITE(6,*) 'EXITING EIRENE TO CHECK TRAPPED PARTICLE'
        CALL EXIT
      ENDIF


c slmod end
      IF (ITIME.EQ.1) THEN
        IF (NLRAD) THEN
          CALL TIMER(TS)
C
          IF (TL.LT.TS.OR.TT.LT.TS) THEN
            MRSURF=0
            IPOLGN=0
C  COLLISION WITH ADDITIONAL SURFACE
            IF (TL.LE.TT) THEN
c slmod begin - debug - tr
              IF (printopt.GE.1.AND.printopt.LE.10)
     .          WRITE(6,'(12X,A)') 'Collision with additional surface'
c slmod end
              ZDT1=TL-ZT
              TL=ZT+ZDT1
              ZTST=TL
              ISRFCL=1
C  COLLISION WITH TIME SURFACE
            ELSEIF (TT.LT.TL) THEN
c slmod begin - debug - tr
              IF (printopt.GE.1.AND.printopt.LE.10)
     .          WRITE(6,'(12X,A)') 'Collision with time surface'
c slmod end
              ZDT1=TT-ZT
              TT=ZT+ZDT1
              ZTST=TT
              ISRFCL=2
            ENDIF
          ELSE
C  COLLISION WITH RADIAL SURFACE
c slmod begin - debug - tr
            IF (printopt.GE.1.AND.printopt.LE.10)
     .        WRITE(6,'(12X,A)') 'Collision with radial surface'
c slmod end
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
          CALL TIMET (ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C
        IF (NLPOL) THEN
          CALL TIMEP(ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C
c slmod begin - debug - tr
c        IF (nacell.GT.1) THEN
c          IF (nacell.EQ.oldnacell) THEN
c            count1 = count1 + 1
c            count2 = count2 + 1
c            IF (count2.EQ.1000) THEN
c              WRITE(0,*) 'OUTSIDE: ',npanu,nacell,NINT(count1)
c              count2 = 0
c            ENDIF
c          ELSE
c            IF (count1.GT.1000) 
c     .        WRITE(0,*) 'OUTSIDE: ',npanu,NINT(oldnacell),NINT(count1)
c            oldnacell = nacell
c            count1    = 1
c            count2    = 1
c          ENDIF
c        ENDIF


        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(3X,2A)')
     .      'FOLNEUT: Output (NR,NP IR,IP MRSURF ZDT1 ',
     .      'X0,Y0,Z0 X00,Y00,Z00)'
          WRITE(6,'(5X,A,2I4,A,2I4,A,I4,E11.3,2X,2(3F9.2,2X))')
     .      '(',NRCELL,NPCELL,') (',IRCELL,IPCELL,')',MRSURF,ZDT1,
     .      X0,Y0,Z0,X00,Y00,Z00

          IF (gridopt.EQ.1) THEN
            ip1 = PolPos(nrcell,x00,y00,cp,'FOLNEUT 1')

            WRITE(6,'(3X,A,I4,3X,A,I4,2F10.3)')
     .        'FOLNEUT: PolPos = ',ip1,' (NR X00,Y00) ',nrcell,x00,y00
            WRITE(6,'(10X,4F14.8)') cp(1),cp(2),cp(3),cp(4)

            IF (ip1.LT.0) THEN
              WRITE(6,'(3X,A)') 'FOLNEUT:'
              WRITE(6,'(3X,A)') 'FOLNEUT: ERROR'
              WRITE(6,'(3X,A)') 'FOLNEUT:'
            ENDIF
          ENDIF
        ENDIF
c slmod end
        IF (ZDT1.LE.0.D0) GOTO 990
      ENDIF
      IF (ZTST.GE.1.D30) GOTO 990
C
C  LOCAL MEAN FREE PATH
C
      IF (IFPATH.NE.1.OR.NRC.LT.0) THEN
        DO 214 J=1,NCOU
          DO 214 K=1,NSTOR
214         XSTOR2(K,J)=0.
        ZMFP=1.D10
      ELSE
        DO 212 J=1,NCOU
          JJ=J
          NCELL=NRCELL+NUPC(J)*NR1P2+NBLCKA
          IF (ITYP.EQ.1) ZMFP=FPATHA(NCELL,CFLAG)
          IF (ITYP.EQ.2) ZMFP=FPATHM(NCELL,CFLAG)

c...sltmp
c          IF (nacell.GT.0) THEN
c            WRITE(6,'(A,3I6,1P,E12.4,0P,I6)') 
c     .        '    NC NA ITYP ZMFP= ',NCELL,nacell,ityp,zmfp,nstor
c            DO i1 = 1, nstor
c              WRITE(6,*) '    XSTOR(K)',xstor(i1)
c            ENDDO
c          ENDIF

          DO 211 K=1,NSTOR
            XSTOR2(K,J)=XSTOR(K)
211       CONTINUE
C  UPDATE INTEGRAL
          ZINT1=ZINT1+CLPD(J)*ZMFPI
          IF (.NOT.NLPR) THEN
            IF (ZINT1.GE.ZLOG) THEN
              IF (NLPOL) NPCELL=NCOUNP(J)
              IF (NLTOR) NTCELL=NCOUNT(J)
              GO TO 213
            ENDIF
            ZINT2=ZINT1
            ZT=ZT+CLPD(J)
          ELSEIF (JCOL.EQ.0) THEN
            IF (ZINT1.GE.ZLOG) THEN
              JCOL=JJ
              IF (NLPOL) NPCELL=NCOUNP(J)
              IF (NLTOR) NTCELL=NCOUNT(J)
            ELSE
              ZINT2=ZINT1
              ZT=ZT+CLPD(J)
            ENDIF
          ENDIF
212     CONTINUE
C
213     CONTINUE
        NCOU=JJ
      ENDIF
C
C  CHECK FOR EVENT
C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(3X,2A)')
     .    'FOLNEUT: Output (NR,NP IR,IP MRSURF ZDT1 ',
     .    'X0,Y0,Z0 X00,Y00,Z00)'
        WRITE(6,'(5X,A,2I4,A,2I4,A,I4,E11.3,2X,2(3F9.2,2X))')
     .    '(',NRCELL,NPCELL,') (',IRCELL,IPCELL,')',MRSURF,ZDT1,
     .    X0,Y0,Z0,X00,Y00,Z00
        WRITE(6,'(3X,A)') 'FOLNEUT: Checking for an event...'
      ENDIF
c slmod end
      IF (NLPR) GOTO 500
      IF (ZINT1.GE.ZLOG) GO TO 220
C
      ZINT2=ZINT1
      ZT=ZTST
C
C  UPDATE CONTRIBUTION TO VOLUME AVERAGED ESTIMATORS
C
      IF (IUPDTE.EQ.1) THEN
        IF (ITYP.EQ.1) CALL UPDATM (XSTOR2)
        IF (ITYP.EQ.2) CALL UPDMOL (XSTOR2)
      ENDIF
C
C  STOP TRACK ?
C
      IF (ISRFCL.EQ.1) CALL ADDCOL (XLI,YLI,ZLI,SG,*104,*380)
      IF (ISRFCL.EQ.2) CALL TIMCOL (AX(2),         *104,*800)
      IF (ISRFCL.EQ.3) CALL TORCOL (               *104)
C
C  NO, CONTINUE TRACK
C
216   CONTINUE
C
C  NEXT CELL - CHECK FOR ESCAPE OR NON DEFAULT ACTING STANDARD SURFACE
C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A)') 'FOLNEUT: NOTHING'
c slmod end
      IF (LEVGEO.NE.4) THEN
C
c      IF (mrsurf.EQ.-1) STOP 'MRSURF.EQ.-1. STOP'

      ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
c slmod begin - new



c...note: Trouble spot
      IF (ists.NE.0) THEN
        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(3X,A,I4)')
     .      'FOLNEUT: ISTS (BEFORE)= ',
     .      ISTS

c Replace with check for NINCS = -1?  Should I id non-exiting surfaces...?
        IF (.TRUE..AND.gridopt.EQ.1.AND.ists  .NE.0     .AND.
     .                  mrsurf .EQ.5.AND.mrsurf.EQ.nrcell) THEN
c          WRITE(6,*) 'WARNING: Non-default surface over-ride - OFF!'
c          WRITE(0,*) 'WARNING: Non-default surface over-ride - OFF!'
          ists = 0
        ENDIF

        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(15X,A,I4,A,2I4)')
     .      '  (AFTER )= ',ISTS,
     .      '  MRSURF,IPCELL= ',mrsurf,ipcell
      ENDIF
c slmod end
      IF (NLRAD.AND.ISTS.NE.0) THEN
        SG=ISIGN(1,NINCX)
        NLSRFX=.TRUE.
        CALL STDCOL (ISTS,1,SG,*104,*380)
      ENDIF
      ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
      IF (NLTOR.AND.ISTS.NE.0) THEN
        SG=ISIGN(1,NINCZ)
        NLSRFZ=.TRUE.
        CALL STDCOL (ISTS,3,SG,*104,*380)
      ENDIF
      ISTS=INMP2I(IRCELL,MPSURF,ITCELL)
      IF (NLPOL.AND.ISTS.NE.0) THEN
        SG=ISIGN(1,NINCY)
        NLSRFY=.TRUE.
        CALL STDCOL (ISTS,2,SG,*104,*380)
      ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A,F12.4)') 'FOLNEUT: A E0=',e0
c slmod end
C
      ELSEIF (LEVGEO.EQ.4) THEN
      ISTS=ABS(INMTI(IPOLGN,MRSURF))
      IF (NLRAD.AND.ISTS.NE.0) THEN
        SG=ISIGN(1,NINCX)
        NLSRFX=.TRUE.
        CALL STDCOL (ISTS,1,SG,*104,*380)
      ENDIF
      ENDIF
C
C
c slmod begin - tr

      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A,I3)')
     .    'FOLNEUT: Crossing radial surface (NINCX) ',nincx

c...DEV:

c slmod end
      NRCELL=NRCELL+NINCX
      IF (NRCELL.GT.NR1STM) GOTO 990
      IF (NACELL.LT.1.AND.NRCELL.LT.1) GOTO 990
c slmod begin - tr
c...note: The above IF statement is new (NACELL was not around before).
      IF (gridopt.EQ.1) THEN
        CALL FindPoloidalCell(NRCELL,X00,Y00)
        IF (NPCELL.EQ.-1) GOTO 990
      ENDIF
c slmod end

C
C  PARTICLE ON SURFACE MRSURF BELONGING TO 1ST (RADIAL OR X-) GRID
C  OR, FOR SOME REASON, IT HAS STOPPED IN THE MIDDLE OF AN ADDITIONAL CELL
      IF (MRSURF.EQ.0) THEN
C  ADVANCE IN SAME ADDITIONAL CELL, AND CONTINUE TRACK
        X0=X0+VELX*ZT
        Y0=Y0+VELY*ZT
        Z0=Z0+VELZ*ZT
        TIME=TIME+ZT/VEL
        IPOLG=IPOLGN
        IF (NLTRA) THEN
          PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
          X01=X0+RMTOR
        ENDIF
        X00=X0
        Y00=Y0
        Z00=Z0
        Z01=Z0
        GOTO 104
      ENDIF
C
C  CHECK IF WE HAVE ENCOUNTERED A SPLITTING ZONE
C  SPLITTING AND RR NOT READY FOR LEVGEO.GE.4
      IF (LEVGEO.LE.3) THEN
      IF (NLSPLT(MRSURF).AND.NLEVEL.LT.MAXLEV.AND.ICOL.EQ.0) THEN
        CALL SPLTRR(1,MRSURF,NINCX,*210,*700)
      ENDIF
      ENDIF
C
C  SWITCH OFF CONDITIONAL EXP. ESTIMATOR ?
      IF (AX(2).LT.WMINC) THEN
        IF (ICOL.EQ.1) GOTO 512
C  NO COLLISION YET; RESTART AGAIN WITH COND. EXP. ESTIMATOR
C                    IN NEW CELL
        AX(1)=1.
        AX(2)=1.
        JCOL=0
      ENDIF
      GOTO 210
C
C  POINT OF COLLISION  220 -- 240
C
220   CONTINUE
c slmod begin - debug - tr
        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(3X,A)') 'FOLNEUT: COLLISION'
c slmod end
      CLPD(NCOU)=(ZLOG-ZINT2)*ZMFP
      ZTC=ZT+CLPD(NCOU)
      IF (IUPDTE.EQ.1) THEN
        IF (ITYP.EQ.1) CALL UPDATM (XSTOR2)
        IF (ITYP.EQ.2) CALL UPDMOL (XSTOR2)
      ENDIF
c slmod begin - tr
      IF ( OPTZMOTION.EQ.0.OR.
     .    (OPTZMOTION.EQ.2.AND.NACELL.NE.0).OR.
     .    (OPTZMOTION.EQ.3.AND.NACELL.NE.0.AND.
     .     (X0.GT.62.5.OR.Y0.LT.-61.0)).OR.
     .    (OPTZMOTION.EQ.4.AND.(NACELL.EQ.0.OR.
     .     (X0.LT.62.5.AND.Y0.GT.-61.0))).OR.
     .    (OPTZMOTION.EQ.5.AND.NACELL.EQ.0)) Z0=Z0+VELZ*ZTC
      X0=X0+VELX*ZTC
      Y0=Y0+VELY*ZTC
c
c      X0=X0+VELX*ZTC
c      Y0=Y0+VELY*ZTC
c      Z0=Z0+VELZ*ZTC
c slmod end
      TIME=TIME+ZTC/VEL
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=LEARC2(X0,Y0,NRCELL,NPANU,'FOLNEUT 2    ')
      ELSEIF (NLFEM) THEN
        IPOLG=0
      ENDIF
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      MSURF=0
      IF (NLTRA) PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
C
230   CONTINUE
C
C  PRE COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL UPCUSR(WS,1)
      ENDIF
C
C  TEST FOR CORRECT CELL NUMBER AT COLLISION POINT
C  KILL PARTICLE, IF TOO LARGE ROUND OFF ERRORS DURING
C  PARTICLE TRACING
C
      IF (NLTEST) CALL CLLTST(*997)
C
C  SAMPLE FROM COLLISION KERNEL FOR NEUTRAL PARTICLES
C  AT PRESENT: NO SUPPRESSION OF ABSORBTION AT IONIZATION
C  FIND NEW WEIGHT, SPECIES INDEX, VELOCITY AND RETURN
C
      IF (ITYP.EQ.1) THEN
        CALL COLATM(CFLAG,COLTYP)
      ELSEIF (ITYP.EQ.2) THEN
        CALL COLMOL(CFLAG,COLTYP)
      ENDIF
      ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
C
C  POST COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL UPCUSR(WS,2)
      ENDIF
C
      IF (COLTYP.EQ.2.) GOTO 700
C
      GOTO 101
C
C  SIMULATION OF COLLISION EVENT FINISHED
C
C
C  ..............................................................
C  .
C  .  INCIDENT ONTO SURFACE
C  ..............................................................
C
380   CONTINUE
C
C
      PR=AX(2)
      IF (ILIIN(MSURF).LE.-2) PR=PR*SG
C
C  UPDATE EFFLUXES ONTO SURFACE AND REFLECT PARTICLE
C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A,F12.4)') 'FOLNEUT: B E0=',e0

      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A)') 'FOLNEUT: Calling ESCAPE'
c slmod end
      CALL ESCAPE (PR,SG,*101,*104,*512)
C
C   GOTO 101: START NEW TRACK OF NEUTRAL PARTICLE
C   GOTO 104: CONTINUE THIS TRACK, TRANSPARENT SURFACE HAS BEEN HIT
C   GOTO 512: RESTORE PREVIOUS COLLISION DATA,
C             CONDITIONAL EXPECTATION ESTIMATOR WAS USED
C
      RETURN
C
C
C  ...................................................
C  .                                                 .
C  .  CONDITIONAL EXPECTATION ESTIMATOR  500 -- 599  .
C  ...................................................
C
C
C  CHECK FOR 1.ST COLLISION ALONG TRACK
500   IF (ICOL.EQ.0.AND.ZINT1.GE.ZLOG) GOTO 505
C
C  UPDATE COND. EXP. ESTIMATOR AND ADVANCE TO NEXT SURFACE
C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A)') 'FOLNEUT: Advance to next surface'
c slmod end
502   CONTINUE
      DO 503 ICOU=1,NCOU
        AX(1)=AX(2)
C  ZMFPI IS EQUIVALENCED TO XSTOR2(NSTOR,ICOU)
        EX=CLPD(ICOU)*XSTOR2(NSTOR,ICOU)
        IF (EX.GE.1.D-10) THEN
          EXPM=EXP(-EX)
          AX(1)=AX(1)*(1.-EXPM)/EX
        ELSE
          EXPM=1.
        ENDIF
        CLPD(ICOU)=CLPD(ICOU)*AX(1)
        AX(2)=AX(2)*EXPM
503   CONTINUE
C
      ZINT2=ZINT1
      ZT=ZTST
C
      IF (IUPDTE.EQ.1) THEN
        IF (ITYP.EQ.1) CALL UPDATM (XSTOR2)
        IF (ITYP.EQ.2) CALL UPDMOL (XSTOR2)
      ENDIF
C
C  STOP TRACK ???
C
      IF (ISRFCL.EQ.1) CALL ADDCOL (XLI,YLI,ZLI,SG,*104,*380)
      IF (ISRFCL.EQ.2) CALL TIMCOL (AX(2),         *104,*800)
      IF (ISRFCL.EQ.3) CALL TORCOL (               *104)
C
C  NO, CONTINUE TRACK
C
      GOTO 216
C
C   SAVE DATA OF FIRST COLLISION ALONG CONDITIONAL TRACK
505   CONTINUE
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A)') 'FOLNEUT: Save conditional track'
c slmod end
      ZMFP=1./XSTOR2(NSTOR,JCOL)
      ZDT1C=(ZLOG-ZINT2)*ZMFP
      ZTC=ZT+ZDT1C
      X0C=X0+VELX*ZTC
      Y0C=Y0+VELY*ZTC
      Z0C=Z0+VELZ*ZTC
      TIMEC=TIME+ZTC/VEL
      NRCLLC=NRCELL
      NPCLLC=NPCELL
      NTCLLC=NTCELL
      NACLLC=NACELL
      NBLCKC=NBLOCK
      NCELLC=NCELL
      ITIMEC=ITIME
      IFPTHC=IFPATH
      IUPDTC=IUPDTE
      VELXC=VELX
      VELYC=VELY
      VELZC=VELZ
      VELC=VEL
      E0C=E0
      GENRC=GENER
      WEIGHC=WEIGHT
      IF (NLTRA)
     .  PHIC=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0C,(X0C+RMTOR)),PI2A)
      DO 510 ISTORE=1,NSTOR
        XSTORC(ISTORE)=XSTOR2(ISTORE,JCOL)
510   CONTINUE
      IF (NLTRC) THEN
        PSAVE=PHI
        PHI=PHIC
        CALL CHCTRC(X0C,Y0C,Z0C,16,11)
        PHI=PSAVE
      ENDIF
      ICOL=1
      GOTO 502
C
C   RESTORE PRE COLLISION DATA AND SAMPLE FROM COLLISION KERNEL
512   X0=X0C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(3X,A,F12.4)') 'FOLNEUT: E E0=',e0
c slmod end
      Y0=Y0C
      Z0=Z0C
      TIME=TIMEC
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      MSURF=0
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      NRCELL=NRCLLC
      NPCELL=NPCLLC
      NTCELL=NTCLLC
      NACELL=NACLLC
      NBLOCK=NBLCKC
      NBLCKA=NSTRD*(NBLOCK-1)+NACELL
      NCELL=NCELLC
      ITIME=ITIMEC
      IFPATH=IFPTHC
      IUPDTE=IUPDTC
      VELX=VELXC
      VELY=VELYC
      VELZ=VELZC
      VEL=VELC
      E0=E0C
      GENER=GENRC
      WEIGHT=WEIGHC
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=LEARC2(X0,Y0,NRCELL,NPANU,'FOLNEUT 3    ')
      ELSEIF (NLFEM) THEN
        IPOLG=0
      ENDIF
      PHI=PHIC
      DO 520 J=1,NSTOR
        XSTOR(J)=XSTORC(J)
520   CONTINUE
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,12)
      ICOL=0
      LGPART=.TRUE.
      GOTO 230
C
700   CONTINUE
C  REGULAR STOP IN SUBR. FOLNEUT, CONTINUE IN SUBR. MCARLO
      RETURN
C
800   CONTINUE
C  REGULAR STOP IN SUBR. FOLNEUT, STOP HISTORY, CENSUS ARRAY FULL
      IF (ICOL.EQ.1.AND..NOT.LGLAST) GOTO 512
      LGPART=.FALSE.
      WEIGHT=0.
      RETURN
C
990   CONTINUE
      CALL LEER(1)
      CALL MASAGE ('ERROR IN FOLNEUT, ZDT1 OR NRCELL OUT OF RANGE  ')
      CALL MASAGE ('PARTICLE IS KILLED                             ')
      WRITE (6,*) 'NPANU,NRCELL,ZDT1,ZTST,TL,TS '
      WRITE (6,*) NPANU,NRCELL,ZDT1,ZTST,TL,TS
      GOTO 995
C
991   CONTINUE
      CALL LEER(1)
      CALL MASAGE ('ERROR IN FOLNEUT, E0TERM IS EQUAL TO ZERO      ')
      CALL MASAGE ('THOMPSON DISTRIBUTION NOT WRITTEN              ')
      WRITE (6,*) 'NPANU ',NPANU
      GOTO 995
C
995   WRITE (6,*) 'MRSURF,MPSURF,MTSURF,MASURF ',
     .             MRSURF,MPSURF,MTSURF,MASURF
      X0ERR=X0+ZT*VELX
      Y0ERR=Y0+ZT*VELY
      Z0ERR=Z0+ZT*VELZ
      IF (NLTRC) THEN
        CALL CHCTRC(X0ERR,Y0ERR,Z0ERR,16,15)
      ELSE
        WRITE (6,*) 'X0,Y0,Z0,ZT ',X0,Y0,Z0,ZT
        WRITE (6,*) 'VELX,VELY,VELZ ',VELX,VELY,VELZ
        WRITE (6,*) 'X0ERR,Y0ERR,Z0ERR ',X0ERR,Y0ERR,Z0ERR
      ENDIF
      GOTO 999
997   CALL MASAGE ('ERROR IN FOLNEUT, DETECTED IN SUBR. CLLTST    ')
      CALL MASAGE ('PARTICLE IS KILLED                            ')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,15)
      GOTO 999
C
998   WRITE (6,*) 'ERROR IN FOLNEUT, SPECIES INDEX OUT OF RANGE '
      WRITE (6,*) ' NPANU,ITYP,IATM,IMOL ',NPANU,ITYP,IATM,IMOL
      GOTO 999
C
999   CONTINUE
      PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      CALL LEER(1)
c slmod begin - debug - tr
      nlost = nlost + 1

      IF (debugopt.NE.0) THEN
        WRITE(0,'(1X,2A)')
     +    'FOLNEUT: Particle lost (NR,NP IR,IP MR PTRASH ETRASH ',
     .    'NAPNU)'

        WRITE(0,'(3X,A,I2,I4,A,I2,I4,A,I4,1X,2F9.3,I7)')
     +    '(',NRCELL,NPCELL,') (',IRCELL,IPCELL,')',MRSURF,
     +    PTRASH(ISTRA),ETRASH(ISTRA),NPANU

        IF (gridopt.EQ.1) THEN
          ip1 = PolPos(nrcell,x0err,y0err,cp,'FOLNEUT 6')

          WRITE(0,'(3X,A,I4,3X,2F10.3)')
     .      '(PolPos X0ERR,Y0ERR) ',
     .      ip1,x0err,y0err

          WRITE(0,'(10X,4F12.6)') cp(1),cp(2),cp(3),cp(4)

c      FUNCTION LEARC1 (X,Y,Z,IPO,IAN,IEN,LOGX,LOGY,NP,TEXT)
C
C   LOGX=TRUE: PARTICLE IS ON A RADIAL SURFACE
C   LOGY=TRUE: PARTICLE IS ON A POLOIDAL SURFACE
C
C   FIND RADIAL MESHPOINT NUMBER LEARC1,
C   (AND POLYGON INDEX IPO, IF NLPLG)
C   IF .NOT.LOGX AND .NOT.LOGY
C     SEARCH IN RADIAL CELLS IAN AND IEN, I.E.
C     SEARCH BETWEEN (!!!) RADIAL SURFACES IAN AND IEN+1
C     THIS SEARCH COVERS THE WHOLE POLOIDAL RANGE
C   IF LOGX
C     SEARCH ON (!!!) RADIAL SURF. IAN FOR POLOIDAL MESH NUMBER IPO
C   IF LOGY
C     SEARCH ON (!!!) POLOIDAL SURF. IAN FOR RADIAL MESH NUMBER LEARC1

            isl1 = -1
            isl2 = -1

c            isl1 = LEARC1(x0err,y0err,z0err,isl2,nrcell,nrcell,
c     .                    .FALSE.,.FALSE.,
c     .                    NPANU,'TEST        ')

c            WRITE(6,'(3X,A,I5)') 'FOLNEUT: LEARC1 = ',isl2

        ENDIF
      ENDIF
c slmod end
      RETURN
      END
C
      SUBROUTINE SPLTRR(IDIM,MS,NINC,*,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMSPL'
      DIMENSION RPST(NPARTC),IPST(MPARTC)
      EQUIVALENCE (RPST(1),X0),(IPST(1),IPOLG)
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
C  RETURN 1: CONTINUE FLIGHT
C  RETURN 2: STOP FLIGHT, BECAUSE OF RUSSIAN ROULETTE
C
C  SPLITTING AND RUSSIAN ROULETTE SURFACE MS
C
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
      IG=SIGN(1.D0,RNUMB(IADD+MS))
C
C  RUSSIAN ROULETTE?
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   SPLTRR: Splitting and Russian roulette'
c slmod end
      IF(NINC*IG.GT.0) GO TO 340
C
C  NO - SPLIT PARTICLE
      X0S=X0
      Y0S=Y0
      Z0S=Z0
      TIMES=TIME
C
      X0=X0+VELX*ZT
      Y0=Y0+VELY*ZT
      Z0=Z0+VELZ*ZT
      TIME=TIME+ZT/VEL
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
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,7)
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
        CALL CHCTRC(XM,YM,ZM,16,8)
      ENDIF
      IF (NLSTOR) CALL STORE(0)
      LGPART=.FALSE.
      RETURN 2
      END
C
C
      SUBROUTINE COLLIDE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  SAMPLE FROM COLLISION KERNEL C
C
C  INPUT:  COMPRT, COMMON BLOCK, CONTAINING ACTUAL PARTICLE PARAMETERS
C          CFLAG,  FLAG FOR POST COLLISION KINETICS
C  OUTPUT: COMPRT, MODIFIED TO POST COLLISION PARTICLE PARAMETERS
C          COLTYP, FLAG: =1 CONTINUE IN CALLING ROUTINE
C                           (FOLNEUT OR FOLION)
C                        =2 EXIT FROM CALLING ROUTINE
C                           EITHER ABSORBTION, OR
C                           TRANSITION NEUTRAL-->ION (IF CALLED
C                           BY FOLNEUT), OR
C                           TRANSITION ION-->NEUTRAL (IF CALLED
C                           BY FOLION)
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMXS'
      INCLUDE 'COUTAU'
      INCLUDE 'CESTIM'
      INCLUDE 'CZT1'
      INCLUDE 'COMUSR'
      INCLUDE 'CCONA'
      INCLUDE 'CRAND'
      INCLUDE 'CSDVI'
      DIMENSION CFLAG(6,3),DUMT(3),DUMV(3)
      SAVE

c slmod begin - temp - not tr (not sure it is necessary)
c...TEMP 
      COMMON /SLTEMP1/ e0sum,e0num,e0avg
      REAL*8           e0sum,e0num,e0avg

      COMMON /SLTEMP2/ CXCNT
      INTEGER          CXCNT

c      REAL*8 e0temp1,e0temp2
c      DATA   /e0temp1=0.0D0,e0temp2=0.0D0/
c slmod end

C
      ENTRY COLATM(CFLAG,COLTYP)
C
C  INCIDENT SPECIES: IOLD
      VELXO=VELX
      VELYO=VELY
      VELZO=VELZ
      VELO=VEL
      V0_PARBO=VEL*(VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL))
      V0_PARBO=V0_PARBO*AMUA*RMASSA(IATM)
      E0O=E0
      WGHTO=WEIGHT
      IOLD=IATM
      NOLD=IATM

      IF (IMETCL(NCELL) == 0) THEN
        NCLMT = NCLMT+1
        ICLMT(NCLMT) = NCELL
        IMETCL(NCELL) = NCLMT
      END IF
C
C  FIRST DECIDE: ELECTRON IMPACT OR ION IMPACT
C
      ZEP1=RANF_EIRENE( )*SIGTOT
      SIGSUM=0.
C
      IF (ZEP1.LE.SIGEIT) THEN
C
C  ELECTRON IMPACT COLLISION:
C
c slmod begin - tr
c...    If we get here, then an atom is being ionised.  Make a call
c       to the USRTIM routine to see if the ion is inside one of the
c       time record triangles:
        IF (TIMNUM.GT.0) CALL USRTI1 
c slmod end
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,2)
        IF (NLSTOR) CALL STORE(2)
C  FIND TYP OF ELECTR. IMPACT COLLISION PROCESS: IRDS
        DO 240 IAEI=1,NAEIIM(IOLD)
          IRDS=LGAEI(IOLD,IAEI)
          SIGSUM=SIGSUM+SIGVEI(IRDS)
          IF (ZEP1.LE.SIGSUM) GOTO 245
240     CONTINUE
        IRDS=LGAEI(IOLD,NAEII(IOLD))
245     CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST-ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED
C
        PTOT=P2NDS(IRDS)
C       PTOTAL=PTOT+PPLDS(IRDS,0)
C  ABSORBED WEIGHT: WEIABS. UPDATE COLLISION ESTIMATOR FOR EAPL
C       WEIABS=WEIGHT*PPLDS(IRDS,0)
C
C  COLLISION ESTIMATOR FOR EAAT, EAPL AND EAEL
        IF (IESTEI(IRDS,3).NE.0) THEN
          EAAT(NCELL)=EAAT(NCELL)-WEIGHT*E0
          EAPL(NCELL)=EAPL(NCELL)+WEIGHT*ESIGEI(IRDS,4)
          EAEL(NCELL)=EAEL(NCELL)+WEIGHT*ESIGEI(IRDS,0)
        ENDIF
C
C  ABSORBTION (INTO BULK SPECIES) IS SUPPRESSED
        WEIGHT=WEIGHT*PTOT
C
C  ARE THERE TEST PARTICLE SECONDARIES AT ALL?
        IF (WEIGHT.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          RETURN
        ENDIF
C
        CALL VELOEI(NCELL,IRDS,VELXO,VELYO,VELZO,VELO)
        XGENER=0.D0
C
C  UPDATE COLLISION ESTIMATORS CONTRIBUTION TO EAAT;EAML;EAIO
        IF (ITYP.EQ.1) THEN
          IF (IESTEI(IRDS,3).NE.0) THEN
            EAAT(NCELL)=EAAT(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.2) THEN
          IF (IESTEI(IRDS,3).NE.0) THEN
            EAML(NCELL)=EAML(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.3) THEN
          IF (IESTEI(IRDS,3).NE.0) THEN
            EAIO(NCELL)=EAIO(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ENDIF
        RETURN
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT) THEN
C
C  CHARGE EXCHANGE:
C
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,4)
C
C   FIND SPECIES INDEX OF CHARGE EXCHANGING BULK ION
        SIGSUM=SIGEIT
        DO 271 IACX=1,NACXIM(IATM)
          IRCX=LGACX(IATM,IACX,0)
          IPLS=LGACX(IATM,IACX,1)
          SIGSUM=SIGSUM+SIGVCX(IRCX)
          IF (ZEP1.LT.SIGSUM) GOTO 272
271     CONTINUE
        IRCX=LGACX(IATM,NACXI(IATM),0)
        IPLS=LGACX(IATM,NACXI(IATM),1)
272     CONTINUE
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?
        FRSTP=N1STX(IRCX,3)
        SCNDP=N2NDX(IRCX,3)
        IF (SCNDP.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          RETURN
        ENDIF
C
C  NEW SPECIES TYPE, INDEX AND ENERGY
C  SUPPRESSION OF ABSORBTION AT CX
C  I.E., NO RANDOM DECISION BETWEEN BULK AND TEST SECONDARIES
        WEIGHT=WEIGHT*SCNDP
        ZEP3=RANF_EIRENE( )*SCNDP
        IF (ZEP3.LE.FRSTP) THEN
C  FOLLOW FIRST SECONDARY, SPEED FROM BULK POPULATION
          ITYP=N1STX(IRCX,1)
          NFLAG=CFLAG(3,1)
          CALL VELOCX(NCELL,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .                NFLAG,IRCX,DUMT,DUMV)

          IF (ITYP.NE.1) STOP 'MARK: NOT JUST ATOMS'

          GOTO (273,274,275),ITYP
C
273       CONTINUE
C  1ST SECONDARY IS ATOM
          IATM=N1STX(IRCX,2)
          IF (IATM.EQ.IOLD) THEN
            XGENER=XGENER+1.D0
          ELSE
            XGENER=0.D0
          ENDIF
          IF (NGENA(IATM).GT.0.AND.XGENER.GT.NGENA(IATM)) THEN
C  UPDATE GENERATION LIMIT TALLIES
C  USE POST COLLISION WEIGHT, VELOCITY AND ENERGY
C  SHOULD MAKE NO DIFFERENCE ON AVERAGE, IF GENERATION LIMIT IS VALID.
C  IF NOT, ONLY THIS GIVES CORRECT BALANCES.
            PGENA(IATM,NCELL)=PGENA(IATM,NCELL)-WEIGHT
            EGENA(IATM,NCELL)=EGENA(IATM,NCELL)-WEIGHT*E0
            V0_PARB=VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL)
            V0_PARB=V0_PARB*VEL*AMUA*RMASSA(IATM)
            VGENA(IATM,NCELL)=VGENA(IATM,NCELL)-WEIGHT*V0_PARB
            LGPART=.FALSE.
            LMETSP(IATM)=.TRUE.
            IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,14)
            ITYP=4
            COLTYP=2
            RETURN
          ENDIF
C
          E0=CVRSSA(IATM)*VELQ
C  NEXT LINES: COLLISION ESTIMATOR FOR CHARGE EXCHANGE NO. IRCX
C  CONSERVE CHARGE IN EACH COLLISION, NOT ONLY ON AVERAGE
          IF (IESTCX(IRCX,1).NE.0) THEN
C  IATMN ATOM SPECIES AFTER CX
            IATMN=IATM
            PAAT(IOLD,NCELL) =PAAT(IOLD,NCELL)-WGHTO
            PAAT(IATMN,NCELL)=PAAT(IATMN,NCELL)+WEIGHT
            PAPL(IPLS,NCELL) =PAPL(IPLS,NCELL)-WEIGHT
            PAEL(NCELL)      =PAEL(NCELL)-WEIGHT
            LMETSP(IOLD)=.TRUE.
            LMETSP(IATMN)=.TRUE.
            LMETSP(NSPAMI+IPLS)=.TRUE.
            IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN ION SPECIES AFTER CX
              IPLSN=N2NDX(IRCX,2)
              PAPL(IPLSN,NCELL)=PAPL(IPLSN,NCELL)+WGHTO
              PAEL(NCELL)      =PAEL(NCELL)+WGHTO
              LMETSP(NSPAMI+IPLSN)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).NE.4) THEN
              GOTO 999
            ENDIF
          ENDIF
          IF (IESTCX(IRCX,3).NE.0) THEN
            EAAT(NCELL)=EAAT(NCELL)-E0O*WGHTO
            EAAT(NCELL)=EAAT(NCELL)+E0*WEIGHT
            EAPL(NCELL)=EAPL(NCELL)-E0*WEIGHT
            IF (N2NDX(IRCX,1).EQ.4) THEN
              EAPL(NCELL)=EAPL(NCELL)+E0O*WGHTO
            ELSE
              GOTO 999
            ENDIF
          ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO COPV
          IF (IESTCX(IRCX,2).NE.0) THEN
            IF (NCPVI.GE.2*NPLSI) THEN
              V0_PARB=VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL)
              V0_PARB=V0_PARB*VEL*AMUA*RMASSA(IATM)
              VPLASP=VXIN(IPLS,NCELL)*BXIN(NCELL)+
     .               VYIN(IPLS,NCELL)*BYIN(NCELL)+
     .               VZIN(IPLS,NCELL)*BZIN(NCELL)
              SIG=SIGN(1.D0,VPLASP)
              IAD=NPLSI+IPLS
C ASSUME: NEW ATOM MOMENTUM IS EQUAL TO INCIDENT ION MOMENTUM
c...sltmp
              STOP 'STOP: COPV BUSINESS 01'
              COPV(IAD,NCELL)=COPV(IAD,NCELL)-WEIGHT*V0_PARB*SIG
              LMETSP(NSPTOT+NADVI+NALVI+NCLVI+IAD)=.TRUE.
              IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN ION SPECIES AFTER CX
                IPLSN=N2NDX(IRCX,2)
                IAD=NPLSI+IPLSN
C ASSUME: NEW ION MOMENTUM IS EQUAL TO INCIDENT ATOM MOMENTUM
c...sltmp
                STOP 'STOP: COPV BUSINESS 02'
                COPV(IAD,NCELL)=COPV(IAD,NCELL)+WGHTO*V0_PARBO*SIG
                LMETSP(NSPTOT+NADVI+NALVI+NCLVI+IAD)=.TRUE.
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
                GOTO 999
              ENDIF
            ENDIF
          ENDIF
          COLTYP=1
          RETURN
274       CONTINUE
C  1ST SECONDARY IS MOLECULE
          IMOL=N1STX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSM(IMOL)*VELQ
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=1
          RETURN
275       CONTINUE
C  1ST SECONDARY IS TEST ION
          IION=N1STX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSI(IION)*VELQ
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=2
          RETURN
        ELSE
C  FOLLOW 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
          ITYP=N2NDX(IRCX,1)
          GOTO (277,278,279) ITYP
C
277       CONTINUE
          IATM=N2NDX(IRCX,2)
C         XGENER= ?
C
          E0=CVRSSA(IATM)*VELO*VELO
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=1
          RETURN
278       CONTINUE
          IMOL=N2NDX(IRCX,2)
          XGENER= 0.D0
C
          E0=CVRSSM(IMOL)*VELO*VELO
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=1
          RETURN
279       CONTINUE
          IION=N2NDX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSI(IION)*VELO*VELO
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=2
          RETURN
        ENDIF
C
C  ELASTIC COLLISION
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGELT) THEN
C
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,3)
C   FIND SPECIES INDEX OF BULK ION COLLISION PARTNER
        SIGSUM=SIGEIT+SIGCXT
        DO 281 IAEL=1,NAELIM(IATM)
          IREL=LGAEL(IATM,IAEL,0)
          IPLS=LGAEL(IATM,IAEL,1)
          SIGSUM=SIGSUM+SIGVEL(IREL)
          IF (ZEP1.LT.SIGSUM) GOTO 282
281     CONTINUE
        IREL=LGAEL(IATM,NAELI(IATM),0)
        IPLS=LGAEL(IATM,NAELI(IATM),1)
282     CONTINUE
C
C  NEW SPECIES INDEX AND ENERGY
C       WEIGHT=WEIGHT*1.
C  FOLLOW SECONDARY, NEW SPEED FROM SUBROUTINE VELOEL
C       ITYP=1
        NFLAG=CFLAG(5,1)
        CALL VELOEL(NCELL,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .              NFLAG,IREL,RMASSA(IOLD))
C
        IATM=IOLD
        E0=CVRSSA(IATM)*VELQ
C  DO NOT UPDATE BGK TALLIES HERE
        IBGK=NPBGKP(IPLS,1)
        IF (IBGK.NE.0) GOTO 300
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION
C  ASSUME, AS BEFORE, NO CHANGE IN SPECIES/TYP
        IF (IESTEL(IREL,1).NE.0) THEN
          PAAT(IOLD,NCELL) =PAAT(IOLD,NCELL)-WGHTO
          PAAT(IATM,NCELL) =PAAT(IATM,NCELL)+WEIGHT
          LMETSP(IOLD)=.TRUE.
          LMETSP(IATM)=.TRUE.
        ENDIF
        IF (IESTEL(IREL,3).NE.0) THEN
          EDEL=E0O*WGHTO-E0*WEIGHT
          EAAT(NCELL)      =EAAT(NCELL)-EDEL
          EAPL(NCELL)      =EAPL(NCELL)+EDEL
        ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO COPV
        IF (IESTEL(IREL,2).NE.0) THEN
          IF (NCPVI.GE.2*NPLSI) THEN
            V0_PARB=VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL)
            V0_PARB=V0_PARB*VEL*AMUA*RMASSA(IATM)
            VDEL=V0_PARBO*WGHTO-V0_PARB*WEIGHT
            VPLASP=VXIN(IPLS,NCELL)*BXIN(NCELL)+
     .             VYIN(IPLS,NCELL)*BYIN(NCELL)+
     .             VZIN(IPLS,NCELL)*BZIN(NCELL)
            VDEL=VDEL*SIGN(1.D0,VPLASP)
            IAD=NPLSI+IPLS
c...sltmp
            STOP 'STOP: COPV BUSINESS 03'
            COPV(IAD,NCELL)=COPV(IAD,NCELL)+VDEL
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+IAD)=.TRUE.
          ENDIF
        ENDIF
300     CONTINUE
        COLTYP=1
        RETURN
C
C  GENERAL ION IMPACT COLLISION: NOT READY
C
      ELSE
C
C
        WRITE (6,*) 'ERROR IN COLATM '
        CALL EXIT
C
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,3)
        DO 261 IAPI=1,NAPIIM(IATM)
C   FIND INDEX OF THAT ION IMPACT COLLISION
          SIGSUM=SIGEIT
          IRPI=LGAPI(IATM,IAPI,0)
          IPLS=LGAPI(IATM,IAPI,1)
          SIGSUM=SIGSUM+SIGVPI(IRPI)
          IF (ZEP1.LT.SIGSUM) GOTO 262
261     CONTINUE
        IRPI=LGAPI(IATM,NAPII(IATM),0)
        IPLS=LGAPI(IATM,NAPII(IATM),1)
262     CONTINUE
        IF (PIOPI(IRPI,0).GT.0) THEN
          DO IIO=1,NIONI
            P=PIOPI(IRPI,IIO)
            IF (P.GT.0.D0) THEN
              IION=IIO
C             E0=E0
              VEL=RSQDVI(IION)*SQRT(E0)
              ITYP=3
              COLTYP=2
              RETURN
            ENDIF
          ENDDO
        ELSEIF (PPLPI(IRPI,0).GT.0) THEN
C  IONIZED INTO BULK SPECIES. STOP THIS TRACK
C  NO SUPPRESSION IF ABSORBTION  AT II
          DO IPL=1,NPLSI
            P=PPLPI(IRPI,IPL)
            IF (P.GT.0.D0) THEN
              IPLS=IPL
              ITYP=4
              LGPART=.FALSE.
              COLTYP=2
              RETURN
            ENDIF
          ENDDO
        ENDIF
C
      ENDIF
C
      ENTRY COLMOL(CFLAG,COLTYP)
C
C  INCIDENT SPECIES: IOLD
      VELXO=VELX
      VELYO=VELY
      VELZO=VELZ
      VELO=VEL
      V0_PARBO=VEL*(VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL))
      V0_PARBO=V0_PARBO*AMUA*RMASSM(IMOL)
      E0O=E0
      WGHTO=WEIGHT
      IOLD=IMOL
      NOLD=NATMI+IMOL

      IF (IMETCL(NCELL) == 0) THEN
        NCLMT = NCLMT+1
        ICLMT(NCLMT) = NCELL
        IMETCL(NCELL) = NCLMT
      END IF
C
C  FIRST DECIDE: ELECTRON IMPACT OR ION IMPACT
C
      ZEP1=RANF_EIRENE( )*SIGTOT
      SIGSUM=0.
C
      IF (ZEP1.LE.SIGEIT) THEN
C
C  ELECTRON IMPACT COLLISION:
C
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,2)
        IF (NLSTOR) CALL STORE(2)
C  FIND TYP OF ELECTR. IMPACT COLLISION PROCESS: IRDS
        DO 340 IMEI=1,NMDSIM(IOLD)
          IRDS=LGMEI(IOLD,IMEI)
          SIGSUM=SIGSUM+SIGVEI(IRDS)
          IF (ZEP1.LE.SIGSUM) GOTO 345
340     CONTINUE
        IRDS=LGMEI(IOLD,NMDSI(IOLD))
345     CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST-ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED
C
        PTOT=P2NDS(IRDS)
C       PTOTAL=PTOT+PPLDS(IRDS,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLDS(IRDS,0)
C
C  COLLISION ESTIMATOR FOR EMML, EMPL AND EMEL
        IF (IESTEI(IRDS,3).NE.0) THEN
          EMML(NCELL)=EMML(NCELL)-WEIGHT*E0
          EMPL(NCELL)=EMPL(NCELL)+WEIGHT*ESIGEI(IRDS,4)
          EMEL(NCELL)=EMEL(NCELL)+WEIGHT*ESIGEI(IRDS,0)
        ENDIF
C
C  ABSORBTION (INTO BULK SPECIES) IS SUPPRESSED
        WEIGHT=WEIGHT*PTOT
C
C  ARE THERE TEST PARTICLE SECONDARIES AT ALL?
        IF (WEIGHT.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          RETURN
        ENDIF
C
        CALL VELOEI(NCELL,IRDS,VELXO,VELYO,VELZO,VELO)
        XGENER=0.D0
C
C  UPDATE COLLISION ESTIMATORS CONTRIBUTION TO EMAT;EMML;EMIO
        IF (ITYP.EQ.1) THEN
          IF (IESTEI(IRDS,3).NE.0) THEN
            EMAT(NCELL)=EMAT(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.2) THEN
          IF (IESTEI(IRDS,3).NE.0) THEN
            EMML(NCELL)=EMML(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.3) THEN
          IF (IESTEI(IRDS,3).NE.0) THEN
            EMIO(NCELL)=EMIO(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ENDIF
        RETURN
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT) THEN
C
C  CHARGE EXCHANGE:
C
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,4)
C
C   FIND SPECIES INDEX OF CHARGE EXCHANGING BULK ION
        SIGSUM=SIGEIT
        DO 371 IMCX=1,NMCXIM(IMOL)
          IRCX=LGMCX(IMOL,IMCX,0)
          IPLS=LGMCX(IMOL,IMCX,1)
          SIGSUM=SIGSUM+SIGVCX(IRCX)
          IF (ZEP1.LT.SIGSUM) GOTO 372
371     CONTINUE
        IRCX=LGMCX(IMOL,NMCXI(IMOL),0)
        IPLS=LGMCX(IMOL,NMCXI(IMOL),1)
372     CONTINUE
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?
        FRSTP=N1STX(IRCX,3)
        SCNDP=N2NDX(IRCX,3)
        IF (SCNDP.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          RETURN
        ENDIF
C
C  NEW SPECIES TYPE, INDEX AND ENERGY
C  SUPPRESSION OF ABSORBTION AT CX
C  I.E., NO RANDOM DECISION BETWEEN BULK AND TEST SECONDARIES
        WEIGHT=WEIGHT*SCNDP
        ZEP3=RANF_EIRENE( )*SCNDP
        IF (ZEP3.LE.FRSTP) THEN
C  FOLLOW FIRST SECONDARY, SPEED FROM BULK POPULATION
          ITYP=N1STX(IRCX,1)
          NFLAG=CFLAG(3,1)
          CALL VELOCX(NCELL,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .                NFLAG,IRCX,DUMT,DUMV)
          GOTO (392,393,394),ITYP
C
392       CONTINUE
C  1ST SECONDARY IS ATOM
          IATM=N1STX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSA(IATM)*VELQ
C  NEXT LINES: COLLISION ESTIMATOR FOR CHARGE EXCHANGE NO. IRCX
C  CONSERVE CHARGE IN EACH COLLISION, NOT ONLY ON AVERAGE
          IF (IESTCX(IRCX,1).NE.0) THEN
C  IATMN ATOM SPECIES AFTER CX
            IATMN=IATM
            PMML(IOLD,NCELL) =PMML(IOLD,NCELL)-WGHTO
            PMAT(IATMN,NCELL)=PMAT(IATMN,NCELL)+WEIGHT
            PMPL(IPLS,NCELL) =PMPL(IPLS,NCELL)-WEIGHT
            PMEL(NCELL)      =PMEL(NCELL)-WEIGHT
            LMETSP(IOLD)=.TRUE.
            LMETSP(IATMN)=.TRUE.
            LMETSP(NSPAMI+IPLS)=.TRUE.
            IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN ION SPECIES AFTER CX
              IPLSN=N2NDX(IRCX,2)
              PMPL(IPLSN,NCELL)=PMPL(IPLSN,NCELL)+WGHTO
              PMEL(NCELL)      =PMEL(NCELL)+WGHTO
              LMETSP(NSPAMI+IPLSN)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).NE.4) THEN
              GOTO 999
            ENDIF
          ENDIF
          IF (IESTCX(IRCX,3).NE.0) THEN
            EMML(NCELL)=EMML(NCELL)-E0O*WGHTO
            EMAT(NCELL)=EMAT(NCELL)+E0*WEIGHT
            EMPL(NCELL)=EMPL(NCELL)-E0*WEIGHT
            IF (N2NDX(IRCX,1).EQ.4) THEN
              EMPL(NCELL)=EMPL(NCELL)+E0O*WGHTO
            ELSE
              GOTO 999
            ENDIF
          ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO COPV
          IF (IESTCX(IRCX,2).NE.0) THEN
            IF (NCPVI.GE.3*NPLSI) THEN
              V0_PARB=VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL)
              V0_PARB=V0_PARB*VEL*AMUA*RMASSM(IMOL)
              VPLASP=VXIN(IPLS,NCELL)*BXIN(NCELL)+
     .               VYIN(IPLS,NCELL)*BYIN(NCELL)+
     .               VZIN(IPLS,NCELL)*BZIN(NCELL)
              SIG=SIGN(1.D0,VPLASP)
              IAD=2*NPLSI+IPLS
C ASSUME: NEW ATOM MOMENTUM IS EQUAL TO INCIDENT ION MOMENTUM
c...sltmp
              STOP 'STOP: COPV BUSINESS 04'
              COPV(IAD,NCELL)=COPV(IAD,NCELL)-WEIGHT*V0_PARB*SIG
              LMETSP(NSPTOT+NADVI+NALVI+NCLVI+IAD)=.TRUE.
              IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN ION SPECIES AFTER CX
                IPLSN=N2NDX(IRCX,2)
                IAD=2*NPLSI+IPLSN
C ASSUME: NEW ION MOMENTUM IS EQUAL TO INCIDENT ATOM MOMENTUM
c...sltmp
                STOP 'STOP: COPV BUSINESS 05'
                COPV(IAD,NCELL)=COPV(IAD,NCELL)+WGHTO*V0_PARBO*SIG
                LMETSP(NSPTOT+NADVI+NALVI+NCLVI+IAD)=.TRUE.
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
                GOTO 999
              ENDIF
            ENDIF
          ENDIF
          COLTYP=1
          RETURN
393       CONTINUE
C  1ST SECONDARY IS MOLECULE
          IMOL=N1STX(IRCX,2)
          IF (IMOL.EQ.IOLD) THEN
            XGENER=XGENER+1.D0
          ELSE
            XGENER=0.D0
          ENDIF
          IF (NGENM(IMOL).GT.0.AND.XGENER.GT.NGENM(IMOL)) THEN
C  UPDATE GENERATION LIMIT TALLIES
            PGENM(IMOL,NCELL)=PGENM(IMOL,NCELL)-WEIGHT
            EGENM(IMOL,NCELL)=EGENM(IMOL,NCELL)-WEIGHT*E0
            V0_PARB=VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL)
            V0_PARB=V0_PARB*VEL*AMUA*RMASSM(IMOL)
            VGENM(IMOL,NCELL)=VGENM(IMOL,NCELL)-WEIGHT*V0_PARB
            LMETSP(NATMI+IMOL)=.TRUE.
            LGPART=.FALSE.
            IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,14)
            ITYP=4
            COLTYP=2
            RETURN
          ENDIF
C
          E0=CVRSSM(IMOL)*VELQ
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=1
          RETURN
394       CONTINUE
C  1ST SECONDARY IS TEST ION
          IION=N1STX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSI(IION)*VELQ
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=2
          RETURN
        ELSE
C  FOLLOW 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
          ITYP=N2NDX(IRCX,1)
          GOTO (395,396,397) ITYP
C
395       CONTINUE
          IATM=N2NDX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSA(IATM)*VELO*VELO
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=1
          RETURN
396       CONTINUE
          IMOL=N2NDX(IRCX,2)
C         XGENER = ?
C
          E0=CVRSSM(IMOL)*VELO*VELO
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=1
          RETURN
397       CONTINUE
          IION=N2NDX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSI(IION)*VELO*VELO
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=2
          RETURN
        ENDIF
C
C  ELASTIC COLLISION
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGELT) THEN
C
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,3)
C   FIND SPECIES INDEX OF BULK ION COLLISION PARTNER
        SIGSUM=SIGEIT+SIGCXT
        DO 398 IMEL=1,NMELIM(IMOL)
          IREL=LGMEL(IMOL,IMEL,0)
          IPLS=LGMEL(IMOL,IMEL,1)
          SIGSUM=SIGSUM+SIGVEL(IREL)
          IF (ZEP1.LT.SIGSUM) GOTO 399
398     CONTINUE
        IREL=LGMEL(IMOL,NMELI(IMOL),0)
        IPLS=LGMEL(IMOL,NMELI(IMOL),1)
399     CONTINUE
C
C  NEW SPECIES INDEX AND ENERGY
C       WEIGHT=WEIGHT*1.
C  FOLLOW SECONDARY, NEW SPEED FROM SUBROUTINE VELOEL
C       ITYP=2
        NFLAG=CFLAG(5,1)
        CALL VELOEL(NCELL,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .              NFLAG,IREL,RMASSM(IOLD))
C
        IMOL=IOLD
        E0=CVRSSM(IMOL)*VELQ
C  DO NOT UPDATE BGK TALLIES HERE
        IBGK=NPBGKP(IPLS,1)
        IF (IBGK.NE.0) GOTO 400
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION
C  ASSUME, AS BEFORE, NO CHANGE IN SPECIES/TYP
        IF (IESTEL(IREL,1).NE.0) THEN
          PMML(IOLD,NCELL) =PMML(IOLD,NCELL)-WGHTO
          PMML(IMOL,NCELL) =PMML(IMOL,NCELL)+WEIGHT
          LMETSP(NATMI+IOLD)=.TRUE.
          LMETSP(NATMI+IMOL)=.TRUE.
        ENDIF
        IF (IESTEL(IREL,3).NE.0) THEN
          EDEL=E0O*WGHTO-E0*WEIGHT
          EMML(NCELL)      =EMML(NCELL)-EDEL
          EMPL(NCELL)      =EMPL(NCELL)+EDEL
        ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO COPV
        IF (IESTEL(IREL,2).NE.0) THEN
          IF (NCPVI.GE.3*NPLSI) THEN
            V0_PARB=VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL)
            V0_PARB=V0_PARB*VEL*AMUA*RMASSM(IMOL)
            VDEL=V0_PARBO*WGHTO-V0_PARB*WEIGHT
            VPLASP=VXIN(IPLS,NCELL)*BXIN(NCELL)+
     .             VYIN(IPLS,NCELL)*BYIN(NCELL)+
     .             VZIN(IPLS,NCELL)*BZIN(NCELL)
            VDEL=VDEL*SIGN(1.D0,VPLASP)
            IAD=2*NPLSI+IPLS
c...sltmp
c            STOP 'STOP: COPV BUSINESS 06'
c slmod begin - new - tr
            COPV (IAD,NCELL     )=COPV (IAD,NCELL     )+VDEL
            COPV2(IAD,NCELL,IN11)=COPV2(IAD,NCELL,IN11)+VDEL
c
c            WRITE(0,*) 'COPV2-11: ',iad,ncell,copv2(iad,ncell,11)
c
c            COPV(IAD,NCELL)=COPV(IAD,NCELL)+VDEL
c slmod end
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+IAD)=.TRUE.
          ENDIF
        ENDIF
400     CONTINUE
        COLTYP=1
        RETURN
C
C  GENERAL ION IMPACT COLLISION: NOT READY
C
      ELSE
C
        WRITE (6,*) 'ERROR IN COLMOL '
        CALL EXIT
C
      ENDIF
C
      ENTRY COLION(CFLAG,COLTYP)
C
C  INCIDENT SPECIES: IOLD
      VELXO=VELX
      VELYO=VELY
      VELZO=VELZ
      VELO=VEL
      V0_PARBO=VEL*(VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL))
      V0_PARBO=V0_PARBO*AMUA*RMASSI(IION)
      E0O=E0
      WGHTO=WEIGHT
      IOLD=IION
      NOLD=NSPAM+IION

      IF (IMETCL(NCELL) == 0) THEN
        NCLMT = NCLMT+1
        ICLMT(NCLMT) = NCELL
        IMETCL(NCELL) = NCLMT
      END IF
C
C  FIRST DECIDE: ELECTRON IMPACT OR ION IMPACT
C
      ZEP1=RANF_EIRENE( )*SIGTOT
      SIGSUM=0.
C
      IF (ZEP1.LE.SIGEIT) THEN
C
C  ELECTRON IMPACT COLLISION:
C
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,2)
        IF (NLSTOR) CALL STORE(2)
C  FIND TYP OF ELECTR. IMPACT COLLISION PROCESS: IRDS
        DO 440 IIDS=1,NIDSIM(IOLD)
          IRDS=LGIEI(IOLD,IIDS)
          SIGSUM=SIGSUM+SIGVEI(IRDS)
          IF (ZEP1.LE.SIGSUM) GOTO 445
440     CONTINUE
        IRDS=LGIEI(IOLD,NIDSI(IOLD))
445     CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST-ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED
C
        PTOT=P2NDS(IRDS)
C       PTOTAL=PTOT+PPLDS(IRDS,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLDS(IRDS,0)
C
C  COLLISION ESTIMATOR FOR EIIO, EIPL AND EIEL
        IF (IESTEI(IRDS,3).NE.0) THEN
          EIIO(NCELL)=EIIO(NCELL)-WEIGHT*E0
          EIPL(NCELL)=EIPL(NCELL)+WEIGHT*ESIGEI(IRDS,4)
          EIEL(NCELL)=EIEL(NCELL)+WEIGHT*ESIGEI(IRDS,0)
        ENDIF
C
C  ABSORBTION (INTO BULK SPECIES) IS SUPPRESSED
        WEIGHT=WEIGHT*PTOT
C
C  ARE THERE TEST PARTICLE SECONDARIES AT ALL?
        IF (WEIGHT.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          RETURN
        ENDIF
C
        CALL VELOEI(NCELL,IRDS,VELXO,VELYO,VELZO,VELO)
C
C  UPDATE COLLISION ESTIMATORS CONTRIBUTION TO EIAT;EIML;EIIO
        IF (ITYP.EQ.1) THEN
          IF (IESTEI(IRDS,3).NE.0) THEN
            EIAT(NCELL)=EIAT(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ELSEIF (ITYP.EQ.2) THEN
          IF (IESTEI(IRDS,3).NE.0) THEN
            EIML(NCELL)=EIML(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ELSEIF (ITYP.EQ.3) THEN
          IF (IESTEI(IRDS,3).NE.0) THEN
            EIIO(NCELL)=EIIO(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ENDIF
        RETURN
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT) THEN
C
C  CHARGE EXCHANGE:
C
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,4)
C
C   FIND SPECIES INDEX OF CHARGE EXCHANGING BULK ION
        SIGSUM=SIGEIT
        DO 490 IICX=1,NICXIM(IION)
          IRCX=LGICX(IION,IICX,0)
          IPLS=LGICX(IION,IICX,1)
          SIGSUM=SIGSUM+SIGVCX(IRCX)
          IF (ZEP1.LT.SIGSUM) GOTO 491
490     CONTINUE
        IRCX=LGICX(IION,NICXI(IION),0)
        IPLS=LGICX(IION,NICXI(IION),1)
491     CONTINUE
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?
        FRSTP=N1STX(IRCX,3)
        SCNDP=N2NDX(IRCX,3)
        IF (SCNDP.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          RETURN
        ENDIF
C
C  NEW SPECIES TYPE, INDEX AND ENERGY
C  SUPPRESSION OF ABSORBTION AT CX
C  I.E., NO RANDOM DECISION BETWEEN BULK AND TEST SECONDARIES
        WEIGHT=WEIGHT*SCNDP
        ZEP3=RANF_EIRENE( )*SCNDP
        IF (ZEP3.LE.FRSTP) THEN
C  FOLLOW FIRST SECONDARY, SPEED FROM BULK POPULATION
          ITYP=N1STX(IRCX,1)
          NFLAG=CFLAG(3,1)
          CALL VELOCX(NCELL,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .                NFLAG,IRCX,DUMT,DUMV)
          GOTO (492,493,494),ITYP
C
492       CONTINUE
C  1ST SECONDARY IS ATOM
          IATM=N1STX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSA(IATM)*VELQ
C  NEXT LINES: COLLISION ESTIMATOR FOR CHARGE EXCHANGE NO. IRCX
C  CONSERVE CHARGE IN EACH COLLISION, NOT ONLY ON AVERAGE
          IF (IESTCX(IRCX,1).NE.0) THEN
C  IATMN ATOM SPECIES AFTER CX
            IATMN=IATM
            PIIO(IOLD,NCELL) =PIIO(IOLD,NCELL)-WGHTO
            PIAT(IATMN,NCELL)=PIAT(IATMN,NCELL)+WEIGHT
            PIPL(IPLS,NCELL) =PIPL(IPLS,NCELL)-WEIGHT
            PIEL(NCELL)      =PIEL(NCELL)-WEIGHT
            LMETSP(IOLD)=.TRUE.
            LMETSP(IATMN)=.TRUE.
            LMETSP(NSPAMI+IPLS)=.TRUE.
            IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN ION SPECIES AFTER CX
              IPLSN=N2NDX(IRCX,2)
              PIPL(IPLSN,NCELL)=PIPL(IPLSN,NCELL)+WGHTO
              PIEL(NCELL)      =PIEL(NCELL)+WGHTO
              LMETSP(NSPAMI+IPLSN)=.TRUE.
            ELSEIF (N2NDX(IRCX,1).NE.4) THEN
              GOTO 999
            ENDIF
          ENDIF
          IF (IESTCX(IRCX,3).NE.0) THEN
            EIIO(NCELL)=EIIO(NCELL)-E0O*WGHTO
            EIAT(NCELL)=EIAT(NCELL)+E0*WEIGHT
            EIPL(NCELL)=EIPL(NCELL)-E0*WEIGHT
            IF (N2NDX(IRCX,1).EQ.4) THEN
              EIPL(NCELL)=EIPL(NCELL)+E0O*WGHTO
            ELSE
              GOTO 999
            ENDIF
          ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO COPV
          IF (IESTCX(IRCX,2).NE.0) THEN
            IF (NCPVI.GE.4*NPLSI) THEN
              V0_PARB=VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL)
              V0_PARB=V0_PARB*VEL*AMUA*RMASSI(IION)
              VPLASP=VXIN(IPLS,NCELL)*BXIN(NCELL)+
     .               VYIN(IPLS,NCELL)*BYIN(NCELL)+
     .               VZIN(IPLS,NCELL)*BZIN(NCELL)
              SIG=SIGN(1.D0,VPLASP)
              IAD=3*NPLSI+IPLS
C ASSUME: NEW ATOM MOMENTUM IS EQUAL TO INCIDENT ION MOMENTUM
c...sltmp
              STOP 'STOP: COPV BUSINESS 07'
              COPV(IAD,NCELL)=COPV(IAD,NCELL)-WEIGHT*V0_PARB*SIG
              LMETSP(NSPTOT+NADVI+NALVI+NCLVI+IAD)=.TRUE.
              IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN ION SPECIES AFTER CX
                IPLSN=N2NDX(IRCX,2)
                IAD=3*NPLSI+IPLSN
C ASSUME: NEW ION MOMENTUM IS EQUAL TO INCIDENT ATOM MOMENTUM
c...sltmp
                STOP 'STOP: COPV BUSINESS 08'
                COPV(IAD,NCELL)=COPV(IAD,NCELL)+WGHTO*V0_PARBO*SIG
                LMETSP(NSPTOT+NADVI+NALVI+NCLVI+IAD)=.TRUE.
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
                GOTO 999
              ENDIF
            ENDIF
          ENDIF
          COLTYP=2
          RETURN
493       CONTINUE
C  1ST SECONDARY IS MOLECULE
          IMOL=N1STX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSM(IMOL)*VELQ
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=2
          RETURN
494       CONTINUE
C  1ST SECONDARY IS TEST ION
          IION=N1STX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSI(IION)*VELQ
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=1
          RETURN
        ELSE
C  FOLLOW 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
          ITYP=N2NDX(IRCX,1)
          GOTO (495,496,497) ITYP
C
495       CONTINUE
          IATM=N2NDX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSA(IATM)*VELO*VELO
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=2
          RETURN
496       CONTINUE
          IMOL=N2NDX(IRCX,2)
          XGENER=0.D0
C
          E0=CVRSSM(IMOL)*VELO*VELO
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=2
          RETURN
497       CONTINUE
          IION=N2NDX(IRCX,2)
C         XGENER = ?

          E0=CVRSSI(IION)*VELO*VELO
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          COLTYP=1
          RETURN
        ENDIF
C
C  ELASTIC COLLISION
C
      ELSE
C
C
        WRITE (6,*) 'ERROR IN COLION '
        CALL EXIT
C
C
      ENDIF
C
C
999   WRITE (6,*) 'ERROR IN COLLIDE '
      WRITE (6,*) 'ITYP ',ITYP,IATM,IMOL,IION,IPLS
      CALL EXIT
      END
C
      SUBROUTINE FOLION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     CHARGED PARTICLE, LAUNCHED AT X0,Y0,Z0 IN CELL NRCELL, IPOLG,
C     NPCELL, NTCELL, NACELL, NBLOCK, IS FOLLOWED
C
C  ON INPUT
C     ITYP=3
C   OUTPUT: ITYP=1  NEXT GENERATION ATOM IATM IS GENERATED
C           ITYP=2  NEXT GENERATION MOLECULE IMOL IS GENERATED
C           ITYP=4  NO NEXT GENERATION PARTICLE IS GENERATED
C                   (PARTICLE ABSORBED IN BULK ION SPECIES)
C           ITYP=0  NO NEXT GENERATION PARTICLE IS GENERATED
C                   (PARTICLE ABSORBED ELSEWHERE)
C
C  DIFFERENCES FROM SUBR. FOLNEUT:
C    1) MOTION ALONG B (VELPAR,....)
C    2) ADDITIONALLY: "FOKKER PLANCK COLLISIONS", ISRFCL=4
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMUSR'
      INCLUDE 'CSPEZ'
      INCLUDE 'COMXS'
      INCLUDE 'CTRIG'
      INCLUDE 'CGRID'
      INCLUDE 'CLGIN'
      INCLUDE 'CCONA'
      INCLUDE 'CZT1'
      INCLUDE 'CUPD'
      INCLUDE 'COUTAU'
      INCLUDE 'CESTIM'
      INCLUDE 'CLOGAU'
      INCLUDE 'COMNNL'
      INCLUDE 'CRAND'
      INCLUDE 'CFPLK'
      DIMENSION CFLAG(6,3),DUMT(3),DUMV(3)
      DIMENSION AX(2),XSTOR(NSTOR)
      DIMENSION XSTOR2(NSTOR,N2ND+N3RD)
      EQUIVALENCE (XSTOR(1),SIGVCX(1))
C
C  NO CONDITIONAL EXPECTATION ESTIMATORS FOR TEST IONS
C
C  ENERGY LOSS FREQUENCY (LANGER APPROXIMATION) (1/SEC)
      TAUEQI(XNI,TI)=8.8E-8*XNI*TI**(-1.5)
C
C  ALL CELL INDICES MUST BE KNOWN AT THIS POINT
C  TENTATIVELY ASSUME: A NEXT GENERATION PARTICLE WILL BE BORN
C
100   LGPART=.TRUE.
      XGENER=0
C
C
C  FIND DIRECTION PARALLEL TO B-FIELD
      NUPC(1)=NPCELL-1+(NTCELL-1)*NP2T3
      NCELL=NRCELL+NUPC(1)*NR1P2+NBLCKA
      IF (NCELL.GT.NSBOX.OR.NCELL.LT.1) GOTO 991
      VCOS=VELX*BXIN(NCELL)+VELY*BYIN(NCELL)+VELZ*BZIN(NCELL)
      SIG=SIGN(1.D0,VCOS)
      VELXS=VELX
      VELYS=VELY
      VELZS=VELZ
      VELS=VEL
C  SET ION ENERGY = PARALLEL ENERGY OF THE IONIZED TEST PARTICLE
C  USE B-FIELD LINE AS TRAJECTORY
      IF (NFOLI(IION).EQ.1) CALL RETUSR(SIG)
      VLXPAR=SIG*BXIN(NCELL)
      VLYPAR=SIG*BYIN(NCELL)
      VLZPAR=SIG*BZIN(NCELL)
      VELPAR=ABS(VEL*VCOS)+EPS12
      E0PAR=CVRSSI(IION)*VELPAR*VELPAR
C
C  FOLLOW MOTION OF TEST IONS ALONG B LINES?
      IF (NFOLI(IION).EQ.-1) THEN
C
C  NO, SIMULATE NEXT INELASTIC COLLISION INSTANTANEOUSLY
C
        IF (IION.LE.0.OR.IION.GT.NIONI) GOTO 998
C  WEIGHT TOO SMALL? STOP HISTORY
        IF (WEIGHT.LT.EPS30) THEN
          LGPART=.FALSE.
          RETURN
        ENDIF
        IF (IFPATH.NE.1) GOTO 993
        LOGION(IION,ISTRA)=.TRUE.
        NCOU=1
        ZMFP=FPATHI(NCELL,CFLAG)
        DO 111 K=1,NSTOR
          XSTOR2(K,1)=XSTOR(K)
111     CONTINUE
        CLPD(1)=ZMFP
        IF (IUPDTE.EQ.1) THEN
          CALL UPDION (XSTOR2)
        ENDIF
        ZTC=0.
C  CARRY OUT INELASTIC COLLISION EVENT, DIRECTLY AT PLACE OF BIRTH
        GOTO 230
      ENDIF
C
C  EACH TEST ION TRACK STARTS AT THIS POINT
C
101   CONTINUE
C     IF (ITYP.EQ.3) THEN
        IF (IION.LE.0.OR.IION.GT.NIONI) GOTO 998
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
      NJUMP=0
      VELX=VELXS
      VELY=VELYS
      VELZ=VELZS
      DO 102 J=1,NR1ST
102     TIMINT(J)=0.0
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
      IF (NLIMII(NCELL).LE.NLIMIE(NCELL)) THEN
        CALL TIMEA1 (MSURF,NCELL,NTCELL,X0,Y0,Z0,TIME,
     .               VLXPAR,VLYPAR,VLZPAR,VELPAR,
c slmod begin - tr
     .               MASURF,XLI,YLI,ZLI,SG,TL,NLTRC,NACELL)
c
c     .               MASURF,XLI,YLI,ZLI,SG,TL,NLTRC)
c slmod end
C       NLPR= :NOT AVAILABLE FOR TEST IONS
        ZTST=TL
        ZDT1=TL
        CLPD(1)=ZDT1
        IF (MASURF.NE.0) ISRFCL=1
      ENDIF
C
C TT: DISTANCE UNTIL NEXT TIMESTEP LIMIT IS REACHED
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
C
C TF: DISTANCE UNTIL NEXT FOKKER PLANCK COLLISION
      TAUI=0.
      DI=0.
      TIDI=0.
      DO IPL=1,NPLSI
        DI=DI+DIIN(IPL,NCELL)
        TIDI=TIDI+DIIN(IPL,NCELL)*TIIN(IPL,NCELL)
        TI=TIDI/(DI+EPS60)
        TAUI=TAUI+TAUEQI(DIIN(IPL,NCELL),TIIN(IPL,NCELL))
      ENDDO
      TAUE=1./TAUI
cdr  stepsize=0.1* vel*taue, i.e. 10 fokker planck collisions per e-time
      TF=TAUE*VELPAR*0.1
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
C  USE PARALLEL VELOCITY, I.E. COMPUTE PARALLEL DISTANCES IN GRID
C  THUS ZT,TS,ZTST,ZDT1,CLPD ETC. ARE PARALLEL DISTANCES
C
      IF (ITIME.EQ.1) THEN
        IF (NLRAD) THEN
          VELXS=VELX
          VELYS=VELY
          VELZS=VELZ
          VELX=VLXPAR
          VELY=VLYPAR
          VELZ=VLZPAR
          CALL TIMER(TS)
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
          CALL TIMET (ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C
        IF (NLPOL) THEN
          CALL TIMEP(ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C
        IF (ZDT1.LE.0.D0) GOTO 990
        VELX=VELXS
        VELY=VELYS
        VELZ=VELZS
      ENDIF
      IF (ZTST.GE.1.D30) GOTO 990
C
C  LOCAL MEAN FREE PATH
C  USE PARALLEL VELOCITY, I.E., COMPUTE PARALLEL MFP
C  BECAUSE CLPD IS THE PARALLEL DISTANCE IN EACH CELL (EXCLUD. GYRO)
C  ETC.. E.G LAMBDA(PARALLEL)=VEL(PARALLEL)/SIGV.
C  THE COLLISION FREQUENCY SIGV, HOWEVER, MUST BE COMPUTED USING THE
C  FULL TEST ION VELOCITY VECTOR, BECAUSE IT MAY DEPEND ON THE RELATIV
C  INTERACTION ENERGY: TO BE WRITTEN
C  FOR INTERACTIONS WITH ELECTRONS THIS IS IRRELEVANT
C
      IF (IFPATH.NE.1.OR.NRC.LT.0) THEN
        DO 214 J=1,NCOU
          DO 214 K=1,NSTOR
214         XSTOR2(K,J)=0.
        ZMFP=1.D10
      ELSE
        VELS=VEL
        VEL=VELPAR
        DO 212 J=1,NCOU
          JJ=J
          NCELL=NRCELL+NUPC(J)*NR1P2+NBLCKA
          ZMFP=FPATHI(NCELL,CFLAG)
          DO 211 K=1,NSTOR
            XSTOR2(K,J)=XSTOR(K)
211       CONTINUE
C  UPDATE INTEGRAL
          ZINT1=ZINT1+CLPD(J)*ZMFPI
C         IF (.NOT.NLPR) THEN
            IF (ZINT1.GE.ZLOG) THEN
              IF (NLPOL) NPCELL=NCOUNP(J)
              IF (NLTOR) NTCELL=NCOUNT(J)
              VEL=VELS
              GOTO 213
            ENDIF
            ZINT2=ZINT1
            ZT=ZT+CLPD(J)
C         ELSEIF (JCOL.EQ.0) THEN
C   CONDITIONAL EXPECTATION ESTIMATOR FOR TEST IONS: TO BE WRITTEN
C         ENDIF
212     CONTINUE
        VEL=VELS
C
213     CONTINUE
        NCOU=JJ
      ENDIF
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
      IF (IUPDTE.EQ.1) THEN
        CALL UPDION(XSTOR2)
      ENDIF
C
C  STOP TRACK ?
C
      IF (ISRFCL.EQ.1) CALL ADDCOL(XLI,YLI,ZLI,SG,*104,*380)
      IF (ISRFCL.EQ.2) CALL TIMCOL(AX(2),         *104,*800)
      IF (ISRFCL.EQ.3) CALL TORCOL(               *104)
      IF (ISRFCL.EQ.4) CALL FPKCOL(               *104,*100)
C
C  NO, CONTINUE TRACK
C
216   CONTINUE
C
C  NEXT CELL - CHECK FOR ESCAPE OR NON DEFAULT ACTING STANDARD SURFACE
C
      IF (LEVGEO.NE.4) THEN
C
      ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
      IF (NLRAD.AND.ISTS.NE.0) THEN
        SG=ISIGN(1,NINCX)
        NLSRFX=.TRUE.
        VELXS=VELX
        VELYS=VELY
        VELZS=VELZ
        VELX=VLXPAR
        VELY=VLYPAR
        VELZ=VLZPAR
        CALL STDCOL (ISTS,1,SG,*104,*380)
      ENDIF
      ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
      IF (NLTOR.AND.ISTS.NE.0) THEN
        SG=ISIGN(1,NINCZ)
        NLSRFZ=.TRUE.
        VELXS=VELX
        VELYS=VELY
        VELZS=VELZ
        VELX=VLXPAR
        VELY=VLYPAR
        VELZ=VLZPAR
        CALL STDCOL (ISTS,3,SG,*104,*380)
      ENDIF
      ISTS=INMP2I(IRCELL,MPSURF,ITCELL)
      IF (NLPOL.AND.ISTS.NE.0) THEN
        SG=ISIGN(1,NINCY)
        NLSRFY=.TRUE.
        VELXS=VELX
        VELYS=VELY
        VELZS=VELZ
        VELX=VLXPAR
        VELY=VLYPAR
        VELZ=VLZPAR
        CALL STDCOL (ISTS,2,SG,*104,*380)
      ENDIF
C
      ELSEIF (LEVGEO.EQ.4) THEN
      ISTS=ABS(INMTI(IPOLGN,MRSURF))
      IF (NLRAD.AND.ISTS.NE.0) THEN
        SG=ISIGN(1,NINCX)
        NLSRFX=.TRUE.
        VELXS=VELX
        VELYS=VELY
        VELZS=VELZ
        VELX=VLXPAR
        VELY=VLYPAR
        VELZ=VLZPAR
        CALL STDCOL (ISTS,1,SG,*104,*380)
      ENDIF
      ENDIF
C
C
c slmod begin - not tr
c      IF (output)
c     .WRITE(0,*) 'MARK: ION TRANSPORT  NINCX =',NINCX
c slmod end
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
CDR
      GOTO 210
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
      IF (IUPDTE.EQ.1) THEN
        CALL UPDION (XSTOR2)
      ENDIF
      X0=X0+VLXPAR*ZTC
      Y0=Y0+VLYPAR*ZTC
      Z0=Z0+VLZPAR*ZTC
      TIME=TIME+ZTC/VELPAR
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=LEARC2(X0,Y0,NRCELL,NPANU,'FOLION 2     ')
      ELSEIF (NLFEM) THEN
        IPOLG=0
      ENDIF
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      MSURF=0
      IF (NLTRA) PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
C
230   CONTINUE
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL UPCUSR(WS,1)
      ENDIF
C
C
C  TEST FOR CORRECT CELL NUMBER AT COLLISION POINT
C  KILL PARTICLE, IF TOO LARGE ROUND OFF ERRORS DURING
C  PARTICLE TRACING
C
      IF (NLTEST) CALL CLLTST(*997)
C
C  SAMPLE FROM COLLISION KERNEL FOR TEST IONS
C  AT PRESENT: NO SUPPRESSION OF ABSORBTION AT IONIZATION
C  FIND NEW WEIGHT, SPECIES INDEX, VELOCITY AND RETURN
C
      CALL COLION(CFLAG,COLTYP)
      ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
C
C  POST COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL UPCUSR(WS,2)
      ENDIF
C
      IF (COLTYP.EQ.2.) GOTO 700
C
      GOTO 100
C
C  SIMULATION OF COLLISION EVENT FINISHED
C
C
C   INCIDENT ONTO SURFACE
C
380   CONTINUE
C
C   REFLECTION FROM  SURFACE
C
      VELX=VELXS
      VELY=VELYS
      VELZ=VELZS
C
C  UPDATE EFFLUXES ONTO SURFACE AND REFLECT PARTICLE
      PR=1.
      IF (ILIIN(MSURF).LE.-2) PR=SG
C
C  FOR NONTRANSPARENT SURFACES:
C  ACCELERATION IN SHEATH IS DONE IN SUBR. ESCAPE
C
      CALL ESCAPE(PR,SG,*100,*104,*996)
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
      CALL LEER(1)
      CALL MASAGE ('ERROR IN FOLION,  ZDT1 OR NRCELL OUT OF RANGE  ')
      CALL MASAGE ('PARTICLE IS KILLED                            ')
      WRITE (6,*) 'NPANU,NRCELL,ZDT1,ZTST ',NPANU,NRCELL,ZDT1,ZTST
      WRITE (6,*) 'TL,TS,ZINT1,ZLOG ',TL,TS,ZINT1,ZLOG
      GOTO 995
991   CONTINUE
      CALL LEER(1)
      CALL MASAGE ('ERROR IN FOLION,  NCELL OUT OF RANGE            ')
      CALL MASAGE ('PARTICLE IS KILLED                            ')
      WRITE (6,*) 'NPANU,NCELL,NRCELL,NPCELL,NTCELL '
      WRITE (6,*)  NPANU,NCELL,NRCELL,NPCELL,NTCELL
      GOTO 995
993   CALL MASAGE ('ERROR IN FOLION,  NO PARTICLE TRACING BUT     ')
      CALL MASAGE ('IFPATH.NE.1. PARTICLE IS KILLED               ')
      WRITE (6,*) 'IION ',IION
      GOTO 999
C
995   WRITE (6,*) 'MRSURF,MPSURF,MTSURF,MASURF ',
     .             MRSURF,MPSURF,MTSURF,MASURF
      X0ERR=X0+ZT*VELX
      Y0ERR=Y0+ZT*VELY
      Z0ERR=Z0+ZT*VELZ
      IF (NLTRC) THEN
        CALL CHCTRC(X0ERR,Y0ERR,Z0ERR,16,15)
      ELSE
        WRITE (6,*) 'X0,Y0,Z0,ZT ',X0,Y0,Z0,ZT
        WRITE (6,*) 'VELX,VELY,VELZ ',VELX,VELY,VELZ
        WRITE (6,*) 'X0ERR,Y0ERR,Z0ERR ',X0ERR,Y0ERR,Z0ERR
      ENDIF
      GOTO 999
996   CALL MASAGE ('ERROR IN FOLION, COND. EXP. ESTIM. NOT IN USE ')
      GOTO 999
997   CALL MASAGE ('ERROR IN FOLION,  DETECTED IN SUBR. CLLTST    ')
      CALL MASAGE ('PARTICLE IS KILLED                            ')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,15)
      GOTO 999
998   WRITE (6,*) 'ERROR IN FOLION, SPECIES INDEX OUT OF RANGE '
      WRITE (6,*) ' NPANU,IION ',NPANU,IION
      GOTO 999
C
999   PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      WEIGHT=0.
      CALL LEER(1)
      RETURN
      END
C
C
      SUBROUTINE CLLTST(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CCONA'
      INCLUDE 'CLOGAU'
      INCLUDE 'CGRID'
      INCLUDE 'CUPD'
C
      IF (NACELL.GT.0) RETURN
C
C  TEST FOR RADIAL CELL INDICES NRCELL, IPOLG
C
      IF (LEVGEO.LE.4) THEN
C
        NTEST0=LEARC1(X0,Y0,Z0,IPOLGT,1,NR1STM,.FALSE.,.FALSE.,NPANU,
     .                'CLLTST      ')
C
      ELSEIF (LEVGEO.EQ.5) THEN
C
        NTEST0=LEAUSR(X0,Y0,Z0)
C
      ENDIF
C
      IF(NTEST0.NE.NRCELL.AND..NOT.NLSRFX) THEN
        WRITE (6,*) 'WRONG CELL-NUMBER IN RADIAL DIRECTION'
        CALL MASJ3 ('NRCELL,NTEST0,NPANU      ',NRCELL,NTEST0,NPANU)
        RETURN 1
      ENDIF
C
C  TEST FOR POLOIDAL CELL INDEX
C
      IF (NLPOL) THEN
        IF (LEVGEO.EQ.1) THEN
          NTEST1=LEARCA(Y0,PSURF,1,NP2ND,1,'CLLTST    ')
        ELSEIF (LEVGEO.EQ.2) THEN
          IF (NLCRC) THEN
            WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
            NTEST1=LEARCA(WINK,PSURF,1,NP2ND,1,'CLLTST    ')
          ELSE
            NTEST1=LEARC2(X0,Y0,NRCELL,NPANU,'CLLTST  ')
          ENDIF
        ELSEIF (LEVGEO.EQ.3) THEN
          NTEST1=IPOLGT
        ELSE
          WRITE (6,*) 'ERROR EXIT IN CLLTST, NLPOL ',NPCELL
          CALL EXIT
        ENDIF
        IF (NTEST1.NE.NPCELL.AND..NOT.NLSRFY) THEN
          WRITE (6,*) 'WRONG CELL-NUMBER IN POLOIDAL DIRECTION'
          CALL MASJ3 ('NPCELL,NTEST1,NPANU     ',NPCELL,NTEST1,NPANU)
          RETURN 1
        ENDIF
      ENDIF
C
C  TEST FOR 3RD CELL INDEX
C
      IF (NLTOR) THEN
C  IN CYLINDER APPROXIMATION
        IF (NLTRZ) THEN
          NTEST2=LEARCA(Z0,ZSURF,1,NTTRA,1,'CLLTST      ')
C  IN TOROIDAL APPROXIMATION
        ELSEIF (NLTRA) THEN
          PHITEST=ATAN2(Z0,(RMTOR+X0))
          PHITEST=PHITEST+ZZONE(NTCELL)
          IF (ABS(PHI-PHITEST).GT.EPS10) THEN
            WRITE (6,*) 'PHI,PHITEST  ',PHI,PHITEST,PHI-PHITEST
            CALL MASJ1 ('NPANU   ',NPANU)
          ENDIF
          NTEST2=LEARCA(PHI,ZSURF,1,NTTRA,1,'CLLTST      ')
C
        ENDIF
        IF (NTEST2.NE.NTCELL.AND..NOT.NLSRFZ) THEN
          WRITE (6,*) 'WRONG CELL-NUMBER IN TOROIDAL DIRECTION'
          CALL MASJ3 ('NTCELL,NTEST2,NPANU     ',NTCELL,NTEST2,NPANU)
          IF (NLTRA) THEN
            ZTESTO=ZFULL*NTCELL
            ZTESTU=ZFULL*(NTCELL-1)
            CALL MASR3 ('ZTESTU,PHI,ZTESTO=   ',ZTESTU,PHI,ZTESTO)
            RETURN 1
          ENDIF
        ENDIF
      ELSEIF (NLTRA) THEN
        PHITEST=ATAN2(Z0,(RMTOR+X0))
        NT=LEARCA(PHI,ZSURF,1,NTTRA,1,'CLLTST      ')
        PHITEST=PHITEST+ZZONE(NT)
        IF (ABS(PHI-PHITEST).GT.EPS10) THEN
          WRITE (6,*) 'PHI,PHITEST  ',PHI,PHITEST,PHI-PHITEST
          CALL MASJ1 ('NPANU   ',NPANU)
        ENDIF
C
        X0T=X0+RMTOR
        Z0T=Z0
        TTT=Z0T/(X0T*TANAL)
        IF (ABS(TTT).GT.1.+EPS10) THEN
          WRITE (6,*) 'WRONG CO-ORDINATES IN TOROIDAL DIRECTION'
          CALL MASJ1 ('NPANU   ',NPANU)
          CALL MASR3 ('X0,Z0,TTT                ',X0,Z0,TTT)
          RETURN 1
        ENDIF
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE STDCOL (ISTS,IDIM,SG,*,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  SG:    =SIGN OF COSINE OF ANGLE OF INCIDENCE
C  ISTS:  =SURFACE INDEX IN NSTSI ARRAYS
C  IDIM:  =INDEX (1,2,3) FOR: RADIAL, POLOIDAL OR TOROIDAL SURFACE
C  RETURN 1:  NO SURFACE TALLIES, FLIGHT CONTINUES
C  RETURN 2:  SURFACE TALLIES,
C             THEN ABSORBTION, REFLECTION MODEL OR CONTINUATION OF FLIGHT
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMUSR'
      INCLUDE 'CLOGAU'
      INCLUDE 'CUPD'
      INCLUDE 'CLGIN'
      INCLUDE 'CGRID'
      INCLUDE 'CPOLYG'
      INCLUDE 'CCONA'
      INCLUDE 'CTRIG'
      INCLUDE 'COUTAU'
      INCLUDE 'CESTIM'
C   COLLISION WITH STANDARD SURFACE NO. MSURF=NLIM+ISTS
C   OF THE RADIAL   (OR X-) GRID: IDIM=1
C   OR THE POLOIDAL (OR Y-) GRID: IDIM=2
C   OR THE TOROIDAL (OR Z-) GRID: IDIM=3
C
C  SAVE DATA OF OLD POINT FOR DIAGNOSTICS
c slmod begin - debug - tr
      COMMON /TRASHCOM/ nlost
      INTEGER           nlost

      INTEGER CHKVAC
      LOGICAL LTRANS,CHKTRA

      LTRANS = .FALSE.

      IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10)
     .  WRITE(6,'(4X,A,2I4,E15.8)')
     .    'STDCOL: ISTS,IDIM,SG= ',
     .    ISTS,IDIM,SG
c slmod end
      X0SA=X0
      Y0SA=Y0
      Z0SA=Z0
      MSURFS=MSURF
      NACLLS=NACELL
C  SET NEW POINT ON NON DEFAULT STANDARD SURFACE ISTS
      if (nlfem) then
        MSURF=ISTS
      else
        MSURF=NLIM+ISTS
      endif
c slmod begin - debug - new

c...  Do not recognize this surface if the particle is inside 
c     the standard mesh.  The surface is only required for re-entry of
c     the particle into the standard mesh:

c...  This turns off the surface when the particle is on the grid and the particle 
c     is on the "wrong" surface:
c
c     IPP/08 Krieger - strange: runtime error because radmap(0) is accessed.
c     That should not be the case. Inserted temporary fix but not sure if
c     that is ok.
c
c     Feb/08 - jde - I am guessing that nrcell=0 must be an error but I am not 
c                    that familiar with EIRENE yet - for the time being I will 
c                    put in Karl's work around
c
      IF (ISWICH(5,MSURF).EQ.3.AND.NACELL.EQ.0.AND.
     .    MRSURF.NE.NRCELL+1.AND.
c...     New for secondary PFZ:
c     .    .NOT.(OPTCONMAP.EQ.1.AND.RADMAP(NRCELL).EQ.0)) THEN
     .    .NOT.(OPTCONMAP.EQ.1.AND.RADMAP(max(1,NRCELL)).EQ.0)) THEN
        WRITE(6,*) 'OVER-RIDE ',nrcell,mrsurf
        ISTS=0
        NLSRFX=.FALSE.
        NLSRFY=.FALSE.
        NLSRFZ=.FALSE.
        RETURN
      ENDIF

      IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10)
     .  WRITE(6,'(4X,A,1P,E15.8,0P,3I5,L2,2I6)')
     .    'STDCOL: MOVING PARTICLE  ZT,MSURF,ISTS,NLIM+ISTS,NLFEM= ',
     .    ZT,MSURF,ISTS,NLIM+ISTS,NLFEM,NACELL,IDIM
c slmod end

c slmod begin - tr (modified)
      IF ( OPTZMOTION.EQ.0.OR.
     .    (OPTZMOTION.EQ.2.AND.NACELL.NE.0).OR.
     .    (OPTZMOTION.EQ.3.AND.NACELL.NE.0.AND.
     .     (X0.GT.62.5.OR.Y0.LT.-61.0)).OR. 
     .    (OPTZMOTION.EQ.4.AND.(NACELL.EQ.0.OR.
     .     (X0.LT.62.5.AND.Y0.GT.-61.0)))) Z0=Z0+VELZ*ZT 
      X0=X0+VELX*ZT
      Y0=Y0+VELY*ZT
c
c      X0=X0+VELX*ZT
c      Y0=Y0+VELY*ZT
c      Z0=Z0+VELZ*ZT
c slmod end
      TIME=TIME+ZT/VEL
      ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
      SCOS=SG
      ICOS=SCOS
      IPOLG=IPOLGN
      IF (NLTRA) THEN
        PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
        X01=X0+RMTOR
      ENDIF
      X00=X0
      Y00=Y0
      Z00=Z0
      Z01=Z0
C
      IF (IDIM.EQ.1) THEN
        NLSRFX=.TRUE.
        NLSRFY=.FALSE.
        NLSRFZ=.FALSE.
        IF (ILIIN(MSURF).LE.0) NRCELL=NRCELL+ICOS
C
C  IN NLFEM AND NLGEN OPTIONS: ALL SURFACES ARE IDIM=1 SURFACES
        IF (NLFEM) NRCELL=MRSURF
        IF (NLGEN) NRCELL=MRSURF
C
      ELSEIF (IDIM.EQ.2) THEN
        NLSRFX=.FALSE.
        NLSRFY=.TRUE.
        NLSRFZ=.FALSE.
        IF (ILIIN(MSURF).LE.0) NPCELL=NPCELL+ICOS
C
        IF (NLPLG.AND.ILIIN(MSURF).GT.0.AND.SCOS.GT.0) IPOLG=IPOLG-1
C
      ELSEIF (IDIM.EQ.3) THEN
        NLSRFX=.FALSE.
        NLSRFZ=.FALSE.
        NLSRFZ=.TRUE.
        IF (ILIIN(MSURF).LE.0) NTCELL=NTCELL+ICOS
      ENDIF
C
      IWEI=ILSIDE(MSURF)*ICOS
      IF (IWEI.LT.0) GOTO 300
C
      IF (ILIIN(MSURF).EQ.2) GOTO 400
C
C  OPERATE A SWITCH
C
      IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10)
     .  WRITE(6,*) '   STDCOL: SWITCH',MSURF,ILSWCH(MSURF)
      IF (ILSWCH(MSURF).NE.0) THEN
C
C  TURN ON OR OFF THE STANDARD GRID CALCULATION
        IF (ISWICH(1,MSURF).NE.0) ITIME=ICOS*ISWICH(1,MSURF)
C  TURN ON OR OFF MFP AND REACTION RATES: PARTICLE ENTERS VACUUM
        IF (ISWICH(2,MSURF).NE.0) IFPATH=ICOS*ISWICH(2,MSURF)
C  TURN ON OR OFF VOLUME AVERAGED TALLIES
        IF (ISWICH(3,MSURF).NE.0) IUPDTE=ICOS*ISWICH(3,MSURF)
C  NEW ADD. CELL INDEX NACELL
C  NOTE: STANDARD SURFACES CANNOT SWITCH NACELL, NBLOCK, NRCELL,....
C       IF (ISWICH(4,MSURF).NE.0) THEN
C  ENTRANCE INTO STANDARD MESH, BLOCK ILBLCK
C  OR
C  EXIT FROM STANDARD MESH, INTO NACELL=ILACLL
        IF (ISWICH(5,MSURF).NE.0) THEN
c slmod begin - new 
          IF (NACELL.EQ.0) THEN
C  SET CELL INDEX EQUAL TO ILACLL
            IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10)
     .        WRITE(6,'(3X,A,2I5)')
     .          ' STDCOL: ENTERING ADDITIONAL SURFACE  MSURF,ILACLL= ',
     .          MSURF,ILACLL(MSURF)

c...GERMANY:  Well, this is to keep the indecies from changing if the particle
c    is bouncing off the back of the target.  It is going to be succeeded by
c    some check of NBLOCK to see if the target
c    should be transparent (divertor gaps):
c              WRITE(80,*) '--->',msurf,iliin(msurf),npanu

c            IF (MSURF.EQ.NLIM+8.AND.EIRNSDTOR.GT.0.AND.
c     .          NBLOCK.EQ.2) LTRANS = .TRUE.

                PHIHOLD=PHI
                PHI=ATAN2(Z0,X0)
            IF (EIRNTRANS.GT.0.AND.CHKTRA(Z0,MSURF)) LTRANS=.TRUE.
                PHI=PHIHOLD

            IF (.NOT.LTRANS.AND.ILIIN(MSURF).GT.0) GOTO 90


c            WRITE(6,*) 'MARK VAC: A ENTERING VACUUM       ',
c     .                 ILACLL(MSURF),NPCELL,NRCELL,npanu
c            WRITE(6,*) 'LEAVING STANDARD MESH ',ILACLL(MSURF),X0,Y0           

           
            IF (OPTVAC.GT.0.AND.ILACLL(MSURF).EQ.998) THEN

c  ASSUME THAT THE SURFACE IS TRANSPARENT!

c  SEARCH ADDITIONAL CELL MESH (VACUUM GRID)
              IF (NLTRA) THEN
                IF (NLTOR) THEN
                ELSE
                  PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
c PHI DOES NOT WORK HERE!  NOT SURE WHY... USING ATAN2(Z0,Z0) INSTEAD
                  TOR=(ATAN2(Z0,X0)+DBLE(NTRSEG-1)*PHISEG)*RADDEG
c...              Special case for (almost) full toroidal grid:
                  IF (NTRSEG.EQ.1.AND.TOR.LT.0.0D0) TOR = TOR + 360.0D0
                ENDIF
              ELSE
                TOR=Z0
              ENDIF

              NACELL=CHKVAC(-1,-1,X0,Y0,TOR,1,NPANU)
c              NACELL=CHKVAC(-1,-1,X0,Y0,Z0,1,NPANU)
	   
              IF (NACELL.EQ.-1.OR.NACELL.EQ.-2) THEN
                WRITE(6,*) 'COULD NOT FIND ADDITIONAL CELL WHEN '//
     .                     'EXITING STANDARD GRID'
                OPTVAC = 2
	   
                IDUM1 = CHKVAC(-1,-1,X0,Y0,TOR,1,NPANU)
c                IDUM1 = CHKVAC(-1,-1,X0,Y0,Z0,1,NPANU)
	   
                OPTVAC = 1
	   
                IWEI=-10
                GOTO 300
              ENDIF
	   
              NBLOCK=NBMLTP
              NRCELL=0
              NPCELL=1
              NTCELL=1
              IF (.NOT.NLADD.OR.NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
                IWEI=-10
                GOTO 300
              ENDIF

            ELSEIF (OPTVAC.GT.0.AND.ILACLL(MSURF).EQ.997) THEN
              IF (NLTRA) THEN
                IF (NLTOR) THEN
                ELSE
                  PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
                  TOR=(ATAN2(Z0,X0)+DBLE(NTRSEG-1)*PHISEG)*RADDEG
c...              Special case for (almost) full toroidal grid:
                  IF (NTRSEG.EQ.1.AND.TOR.LT.0.0D0) TOR = TOR + 360.0D0
                ENDIF
              ELSE
                TOR=Z0
              ENDIF
c  SEARCH ADDITIONAL CELL MESH (VACUUM GRID)
              NACELL=CHKVAC(-1,-1,X0,Y0,TOR,2,NPANU)
c              NACELL=CHKVAC(-1,-1,X0,Y0,Z0,2,NPANU)
              IF (NACELL.EQ.-1.OR.NACELL.EQ.-2) THEN
                WRITE(6,*) 'COULD NOT FIND ADDITIONAL CELL WHEN '//
     .                     'EXITING SOL SURFACE OF MAIN GRID, '//
     .                     'ASSUMING 1ST ADDITIONAL CELL'
                NACELL=1
              ENDIF
              NBLOCK=NBMLTP
              NRCELL=0
              NPCELL=1
              NTCELL=1
              IF (.NOT.NLADD.OR.NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
                IWEI=-10
                GOTO 300
              ENDIF
            ELSE
              NACELL=ILACLL(MSURF) 
              NBLOCK=NBMLTP
              NRCELL=0
              NPCELL=1
              NTCELL=1
              IF (.NOT.NLADD.OR.NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
                IWEI=-10
                GOTO 300
              ENDIF
            ENDIF
c
c          IF (NACELL.EQ.0) THEN
cC  SET CELL INDEX EQUAL TO ILACLL
c            NACELL=ILACLL(MSURF)
c            NBLOCK=NBMLTP
c            NRCELL=0
c            NPCELL=1
c            NTCELL=1
c            IF (.NOT.NLADD.OR.NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
c              IWEI=-10
c              GOTO 300
c            ENDIF
c slmod end
c slmod begin - new
          ELSEIF (NACELL.GT.0) THEN
C  ENTRANCE INTO STANDARD MESH, INTO NBLOCK=ILBLCK

            IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10)
     .        WRITE(6,'(3X,A,2I5)')
     .          ' STDCOL: ENTERING STANDARD MESH       MSURF,ILBLCK= ',
     .          MSURF,ILBLCK(MSURF)


c...GERMANY:  Well, this is to keep the indecies from changing if the particle
c    is bouncing off the back of the target.  It is going to be succeeded by
c    some check of NBLOCK to see if the target
c    should be transparent (divertor gaps):
c              WRITE(0,*) '--->',msurf,iliin(msurf),npanu

c            IF (MSURF.EQ.NLIM+8.AND.EIRNSDTOR.GT.0) LTRANS = .TRUE.

                PHIHOLD=PHI
                PHI=ATAN2(Z0,X0)
            IF (EIRNTRANS.GT.0.AND.CHKTRA(Z0,MSURF)) LTRANS=.TRUE.
                PHI=PHIHOLD

            IF (.NOT.LTRANS.AND.ILIIN(MSURF).GT.0) GOTO 90

            IF (NLMLT.AND.ILBLCK(MSURF).EQ.999) THEN
              PHIHOLD=PHI
              PHI=ATAN2(Z0,X0)
              CALL CHKSTD(Z0,NBLOCK)
              PHI=PHIHOLD
            ELSE
              NBLOCK=ILBLCK(MSURF)
            ENDIF
c
c          ELSEIF (NACELL.GT.0) THEN
cC  ENTRANCE INTO STANDARD MESH, INTO NBLOCK=ILBLCK
c            NBLOCK=ILBLCK(MSURF)
c slmod end
            NACELL=0
C  FIND  NRCELL,IPOLG IN STANDARD MESH, BLOCK NBLOCK
            IF (IDIM.EQ.1) THEN
c slmod begin - grid - tr
              IF (GRIDOPT.EQ.1) THEN
c...note: Ambiguous - fix
c               WRITE(6,*) 'NRCELL SHIFT',MRSURF,OPTCONMAP,RADMAP(MRSURF)
                IF (MRSURF.GT.DIVSUR.AND.
     .              .NOT.(OPTCONMAP.EQ.1.AND.RADMAP(MRSURF).EQ.0)) THEN
                   NRCELL = MRSURF - 1
c                  WRITE(6,*) 'USUAL',NRCELL,MRSURF
                ELSE
                  NRCELL=MIN0(NR1STM,MRSURF)
c                  WRITE(6,*) 'NONE',NRCELL,MRSURF
                ENDIF   

c                IF (MRSURF.GT.DIVSUR) NRCELL = MRSURF - 1
              ELSE
                NRCELL=MIN0(NR1STM,MRSURF)
              ENDIF
c
c              NRCELL=MIN0(NR1STM,MRSURF)
c slmod end
              IAN=MRSURF
              NDUM=LEARC1(X0,Y0,Z0,IPOLG,IAN,IEN,NLSRFX,NLSRFY,NPANU,
     .                   'STDCOL      ')
c slmod begin - debug - tr
              IF (printopt.GE.1.AND.printopt.LE.10)
     .          WRITE(6,'(4X,A,3I4)')
     .            'STDCOL: Particle location (MRSURF,NRCELL,IPLOG) ',
     .            mrsurf,nrcell,ipolg
c slmod end
            ELSEIF (IDIM.EQ.2) THEN
              IPOLG=MIN0(NP2NDM,MPSURF)
              IAN=MPSURF
              NRCELL=LEARC1(X0,Y0,Z0,IDUM,IAN,IEN,NLSRFX,NLSRFY,NPANU,
     .                     'STDCOL      ')
            ELSE
              NRCELL=LEARC1(X0,Y0,Z0,IPOLG,1,NR1STM,
     .                     .FALSE.,.FALSE.,NPANU,
     .                     'STDCOL      ')
            ENDIF
C  FIND NTCELL IN STANDARD MESH, BLOCK NBLOCK
            IF (NLTOR) THEN
              IF (NLTRZ) THEN
                NTCELL=LEARCA(Z0,ZSURF,1,NT3RD,1,'STDCOL    ')
              ELSEIF (NLTRA) THEN
                NTCELL=LEARCA(PHI,ZSURF,1,NT3RD,1,'STDCOL    ')
              ENDIF
            ELSE
              NTCELL=1
            ENDIF
C  FIND NPCELL IN STANDARD MESH, BLOCK NBLOCK
            IF (NLPOL) THEN
              IF (LEVGEO.EQ.1) THEN
                NPCELL=LEARCA(Y0,PSURF,1,NP2ND,1,'STDCOL')
              ELSEIF (LEVGEO.EQ.2) THEN
                IF (NLCRC) THEN
                  WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
                  NPCELL=LEARCA(WINK,PSURF,1,NP2ND,1,'STDCOL')
                ELSE
                  NPCELL=LEARC2(X0,Y0,NRCELL,NPANU,'STDCOL')
                ENDIF
              ELSEIF (LEVGEO.EQ.3) THEN
                NPCELL=IPOLG
              ELSE
                WRITE (6,*) 'ERROR EXIT IN STDCOL, NLPOL ',NPCELL
                CALL EXIT
              ENDIF
            ELSE
              NPCELL=1
            ENDIF
c slmod begin - not tr
c            WRITE(6,*) 'MARK VAC: ENTERING STANDARD MESH',
c     .                  NPCELL,NRCELL,npanu
c slmod end
          ENDIF

          NBLCKA=NSTRD*(NBLOCK-1)+NACELL
c slmod begin - tr
90        CONTINUE
c slmod end
C  ENTRANCE INTO STANDARD MESH, BLOCK ILBLCK
C  OR
C  EXIT FROM STANDARD MESH, INTO NACELL=NBLOCK+ILACLL
        ELSEIF (ISWICH(6,MSURF).NE.0) THEN
          IF (NACELL.EQ.0) THEN
C  SET CELL INDEX EQUAL TO ILACLL
            NACELL=NBLOCK+ICOS*ISWICH(6,MSURF)*ILACLL(MSURF)
            NBLOCK=NBMLTP
C
            NRCELL=0
            NPCELL=1
            NTCELL=1
            IF (.NOT.NLADD.OR.NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
              IWEI=-10
              GOTO 300
            ENDIF
          ELSEIF (NACELL.GT.0) THEN
            STOP 'MARK: STOP: NOT ALLOWED TO ENTER STANDARD MESH HERE'
C  ENTRANCE INTO STANDARD MESH, INTO NBLOCK=NACELL+ILBLCK
            NBLOCK=NACELL+ICOS*ISWICH(6,MSURF)*ILBLCK(MSURF)
            NACELL=0
C  FIND  NRCELL,IPOLG IN STANDARD MESH, BLOCK NBLOCK
            IF (IDIM.EQ.1) THEN
              NRCELL=MIN0(NR1STM,MRSURF)
              IAN=MRSURF
              NDUM=LEARC1(X0,Y0,Z0,IPOLG,IAN,IEN,NLSRFX,NLSRFY,NPANU,
     .                   'STDCOL      ')
            ELSEIF (IDIM.EQ.2) THEN
              IPOLG=MIN0(NP2NDM,MPSURF)
              IAN=MPSURF
              NRCELL=LEARC1(X0,Y0,Z0,IDUM,IAN,IEN,NLSRFX,NLSRFY,NPANU,
     .                     'STDCOL      ')
            ELSE
              NRCELL=LEARC1(X0,Y0,Z0,IPOLG,1,NR1STM,
     .                     .FALSE.,.FALSE.,NPANU,
     .                     'STDCOL      ')
            ENDIF
C  FIND NTCELL IN STANDARD MESH, BLOCK NBLOCK
            IF (NLTOR) THEN
              IF (NLTRZ) THEN
                NTCELL=LEARCA(Z0,ZSURF,1,NT3RD,1,'STDCOL    ')
              ELSEIF (NLTRA) THEN
                NTCELL=LEARCA(PHI,ZSURF,1,NT3RD,1,'STDCOL    ')
              ENDIF
            ELSE
              NTCELL=1
            ENDIF
C  FIND NPCELL IN STANDARD MESH, BLOCK NBLOCK
            IF (NLPOL) THEN
              IF (LEVGEO.EQ.1) THEN
                NPCELL=LEARCA(Y0,PSURF,1,NP2ND,1,'STDCOL')
              ELSEIF (LEVGEO.EQ.2) THEN
                IF (NLCRC) THEN
                  WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
                  NPCELL=LEARCA(WINK,PSURF,1,NP2ND,1,'STDCOL')
                ELSE
                  NPCELL=LEARC2(X0,Y0,NRCELL,NPANU,'STDCOL')
                ENDIF
              ELSEIF (LEVGEO.EQ.3) THEN
                NPCELL=IPOLG
              ELSE
                WRITE (6,*) 'ERROR EXIT IN STDCOL, NLPOL ',NPCELL
                CALL EXIT
              ENDIF
            ELSE
              NPCELL=1
            ENDIF
          ENDIF
          NBLCKA=NSTRD*(NBLOCK-1)+NACELL
        ENDIF
      ENDIF
C
C  SWITCHING DONE
C
      IF (NLTRC) THEN
        CALL CHCTRC(X0,Y0,Z0,16,6)
      ENDIF
C
C
c slmod begin - tr 
c... GERMANY:
c...  Surface is only transparent for certain values of NBLOCK:
      IF (LTRANS) RETURN 1
c slmod end
      IF (ILIIN(MSURF).LT.0) THEN
        IF (ILIIN(MSURF).EQ.-1) RETURN 1
        RETURN 2
      ENDIF
C
C  ILIIN(MSURF) .GT. 0, PREPARE REFLECTION, I.E.: SET SURFACE NORMAL
C
      GOTO (100,150,200),IDIM
C
      ENTRY STDNOR(X0E,Y0E,Z0E,IDIM,SCOSE,MSURFE,*,*)
      X0=X0E
      Y0=Y0E
      Z0=Z0E
      SCOS=SCOSE
      MSURF=MSURFE
      IRCELL=NRCELL
      IPCELL=NPCELL
      ITCELL=NTCELL
C
      GOTO (100,150,200,250),IDIM
C
C  RADIAL SURFACE
C
100   CONTINUE
C
      IF (LEVGEO.EQ.3) THEN
        IF (NLPOL)      IP=IPCELL
        IF (.NOT.NLPOL) IP=IPOLG
        CRTX=PLNX(MRSURF,IP)*SCOS
        CRTY=PLNY(MRSURF,IP)*SCOS
        CRTZ=0.
      ELSEIF (LEVGEO.EQ.2) THEN
        CRTX=(X0-EP1(MRSURF))*ELLQ(MRSURF)
        CRTY=Y0
        PHINM=SQRT(CRTX*CRTX+CRTY*CRTY)
        CRTX=CRTX/PHINM*SCOS
        CRTY=CRTY/PHINM*SCOS
        CRTZ=0.
      ELSEIF (LEVGEO.EQ.4) THEN
        IP=IPOLG
        CRTX=PTRIX(IP,MRSURF)*SCOS
        CRTY=PTRIY(IP,MRSURF)*SCOS
        CRTZ=0.
      ELSEIF (LEVGEO.EQ.1) THEN
        CRTX=SCOS
        CRTY=0.
        CRTZ=0.
C  PERIODICITY SURFACE IN X DIRECTION
        IF (ILIIN(MSURF).GT.3) THEN
          M=IDEZ(ILIIN(MSURF),2,2)
C  NEW X0, AND KEEP Y AND Z (NLTRZ) OR Y AND PHI (NLTRA) CONSTANT
          IF (NLTRA) THEN
            TANPHI=Z0/(X0+RMTOR)
            X0=RSURF(M)
            Z0=TANPHI*(X0+RMTOR)
          ELSEIF (NLTRZ) THEN
            X0=RSURF(M)
          ELSEIF (NLTRT) THEN
            WRITE (6,*) 'EXIT IN STDNOR'
            CALL EXIT
          ENDIF
C  NEW CELL NUMBERS
          MRSURF=M
          IF (SCOS.GT.0) NRCELL=M
          IF (SCOS.LT.0) NRCELL=M-1
C
          IF (NLTRC) THEN
            CALL CHCTRC(X0,Y0,Z0,0,9)
          ENDIF
        ENDIF
      ELSEIF (LEVGEO.EQ.5) THEN
C
C  GENERAL GEOMETRY OPTION: PROVIDE OUTER SURFACE NORMAL UNIT VECTOR
C                           CRTX,CRTY,CRTZ
C  IN CASE OF PERIODICITY: ALSO NEW POSITION, SPEED, SURFACE- AND CELL NUMBERS
C
        ISTS=MSURF-NLIM
        CALL NORUSR(ISTS,X0,Y0,Z0,CRTX,CRTY,CRTZ,SCOS,NRCELL)
      ENDIF
      RETURN 2
C
C  POLOIDAL SURFACE
C
150   CONTINUE
C
      IF (LEVGEO.EQ.3) THEN
        IF (NLRAD)      IR=IRCELL
C       IF (.NOT.NLRAD) IR=???
        CRTX=PPLNX(IR,MPSURF)*SCOS
        CRTY=PPLNY(IR,MPSURF)*SCOS
        CRTZ=0.
      ELSEIF (LEVGEO.EQ.2) THEN
        WRITE (6,*) 'STDCOL IDIM=2, LEVGEO=2: TO BE WRITTEN '
        CALL EXIT
      ELSEIF (LEVGEO.EQ.1) THEN
        CRTX=0.
        CRTY=SCOS
        CRTZ=0.
C  PERIODICITY SURFACE IN Y DIRECTION
        IF (ILIIN(MSURF).GT.3) THEN
          M=IDEZ(ILIIN(MSURF),2,2)
          Y0=PSURF(M)
C  NEW CELL NUMBERS
          MPSURF=M
          IF (SCOS.GT.0) NPCELL=M
          IF (SCOS.LT.0) NPCELL=M-1
C
          IF (NLTRC) THEN
            CALL CHCTRC(X0,Y0,Z0,0,9)
          ENDIF
        ENDIF
      ELSE
        WRITE (6,*) 'STDCOL IDIM=2, BUT LEVGEO.GT.3 CALL EXIT '
        CALL EXIT
      ENDIF
      RETURN 2
C
C  TOROIDAL SURFACE
C
200   CONTINUE
C
      IF (NLTRZ) THEN
        CRTX=0.
        CRTY=0.
        CRTZ=SCOS
C  PERIODICITY SURFACE IN Z-DIRECTION
        IF (ILIIN(MSURF).GT.3) THEN
          M=IDEZ(ILIIN(MSURF),2,2)
          Z0=ZSURF(M)
          MTSURF=M
          IF (SCOS.GT.0) NTCELL=M
          IF (SCOS.LT.0) NTCELL=M-1
          IF (NLTRC) THEN
            CALL CHCTRC(X0,Y0,Z0,0,9)
          ENDIF
        ENDIF
      ELSEIF (NLTRA) THEN
        IF (ILIIN(MSURF).GT.3) THEN
          M=IDEZ(ILIIN(MSURF),2,2)
C  NEW POSITION. KEEP X0,Y0 FIXED
          PHI=ZSURF(M)
C  ROTATE VELOCITY: TO BE WRITTEN
C         Z0=?
C         VELX=?
C         VELZ=?
          WRITE (6,*) 'STDCOL IDIM=3, NLTRA, PERIODIC. CALL EXIT '
          CALL EXIT
C  NEW CELL NUMBERS
          MTSURF=M
          IF (SCOS.GT.0) NTCELL=M
          IF (SCOS.LT.0) NTCELL=M-1
C
          IF (NLTRC) THEN
            CALL CHCTRC(X0,Y0,Z0,0,9)
          ENDIF
        ELSE
C         CRTX=?
          CRTY=0.
C         CRTZ=?
          WRITE (6,*) 'STDCOL IDIM=3, BUT NOT PERIODIC. CALL EXIT '
          CALL EXIT
        ENDIF
      ELSE
        WRITE (6,*) 'STDCOL IDIM=3, BUT NEITHER NLTRZ NOR NLTRA: EXIT '
        CALL EXIT
      ENDIF
      RETURN 2
C
C  MIXED SURFACE
C
250   CONTINUE
      IF (NLSRFX) GOTO 100
      IF (NLSRFY) GOTO 150
      WRITE (6,*) ' ERROR IN STDNOR,',
     .            ' IDIM=4 AND NEITHER NLSRFX NOR NLSRFY '
      CALL EXIT
C
C
300   CONTINUE
C
C  IWEI.LT.0, I.E., ILIIN OPTION IS OVERWRITTEN FROM THIS SIDE
C
      IF (IWEI.EQ.-1) THEN
C  PARTICLE HAS HITTEN A SURFACE FROM AN ABSORBING SIDE
C  UPDATE FLUXES (DO NOT SET WEIGHT=0.D0) AND ABSORB PARTICLE
        IF (NLTRC) THEN
          CALL CHCTRC(X0,Y0,Z0,16,6)
          WRITE (6,*) 'ABSORB PARTICLE: NPANU ',NPANU
        ENDIF
        SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        LGPART=.FALSE.
        RETURN 2
      ELSEIF (IWEI.EQ.-2) THEN
C  KILL THIS PARTICLE BECAUSE IT COMES FROM WRONG SIDE
C  DO NOT UPDATE FLUXES (SET WEIGHT=0.D0)
        PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
        ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
        IF (NLTRC) THEN
          CALL CHCTRC(X0,Y0,Z0,16,15)
        ENDIF
        WRITE (6,*) 'ERROR DETECTED IN SUBR. STDCOL '
        WRITE (6,*) 'PARTICLE COMES FROM WRONG SIDE '
        CALL MASJ1 ('NPANU=  ',NPANU)
        CALL MASJ1 ('MSURF NW',-(MSURF-NLIM))
        IF (MSURFS.GT.NLIM) MSURFS=-(MSURFS-NLIM)
        CALL MASJ1 ('MSURF OD',MSURFS)
        CALL MASR3 ('X0,Y0,Z0 (NEW)          ',X0,Y0,Z0)
        CALL MASR3 ('X0,Y0,Z0 (OLD)          ',X0SA,Y0SA,Z0SA)
        CALL MASR3 ('VELX,VELY,VELZ          ',VELX,VELY,VELZ)
        CALL MASR2 ('WEIGHT,E0       ',WEIGHT,E0)
        SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        WEIGHT=0.
        LGPART=.FALSE.
c slmod begin - tr
        NLOST=NLOST+1
c slmod end
        RETURN 2
      ELSEIF (IWEI.EQ.-3) THEN
C  SURFACE IS NOT SEEN BY THE PARTICLE BECAUSE OF ILSIDE OPTION
C  I.E., SURFACE IS TRANSPARENT FROM THIS SIDE
C  ACTS AS ILIIN=0 OPTION (NO SURFACE TALLIES, NO SWITCHES)
        IF (IDIM.EQ.1.AND.ILIIN(MSURF).GT.0) NRCELL=NRCELL+ICOS
        IF (IDIM.EQ.2.AND.ILIIN(MSURF).GT.0) NPCELL=NPCELL+ICOS
        IF (IDIM.EQ.3.AND.ILIIN(MSURF).GT.0) NTCELL=NTCELL+ICOS
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,6)
        RETURN 1
      ELSEIF (IWEI.EQ.-10) THEN
C  KILL THIS PARTICLE BECAUSE CELL NUMBER OUT OF RANGE DUE TO SWITCHING
C  DO NOT UPDATE FLUXES (SET WEIGHT=0.D0)
        PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
        ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
        IF (NLTRC) THEN
          CALL CHCTRC(X0,Y0,Z0,16,15)
        ENDIF
        WRITE (6,*) 'ERROR DETECTED IN SUBR. STDCOL '
        WRITE (6,*) 'NACELL OUT OF RANGE '
        CALL MASJ1 ('NPANU=  ',NPANU)
        WRITE (6,*) 'NLMLT,NLADD ',NLMLT,NLADD
        WRITE (6,*) 'NBMLT,NRADD ',NBMLT,NRADD
        CALL MASJ1 ('MSURF NW',-(MSURF-NLIM))
        IF (MSURFS.GT.NLIM) MSURFS=-(MSURFS-NLIM)
        CALL MASJ1 ('MSURF OD',MSURFS)
        CALL MASJ1 ('NACL NEW',NACELL)
        CALL MASJ1 ('NACL OLD',NACLLS)
        CALL MASR3 ('X0,Y0,Z0 (NEW)          ',X0,Y0,Z0)
        CALL MASR3 ('X0,Y0,Z0 (OLD)          ',X0SA,Y0SA,Z0SA)
        CALL MASR3 ('VELX,VELY,VELZ          ',VELX,VELY,VELZ)
        CALL MASR2 ('WEIGHT,E0       ',WEIGHT,E0)
        SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        WEIGHT=0.
        LGPART=.FALSE.
c slmod begin - tr
        NLOST=NLOST+1
c slmod end
        RETURN 2
      ENDIF
C
C  ABSORBING SURFACE
C  UPDATE FLUXES (DO NOT SET WEIGHT=0.D0)
C
400   CONTINUE
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,6)
      RETURN 2
      END
C
C
      SUBROUTINE FPKCOL(*,*)
C  FOKKER PLANCK ELASTIC COLLISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  RETURN 1:  NOT IN USE
C  RETURN 2:  START COMPLETELY NEW TEST ION TRACK, SAME SPECIES
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMUSR'
      INCLUDE 'COMXS'
      INCLUDE 'CUPD'
      INCLUDE 'CLGIN'
      INCLUDE 'CLOGAU'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
      INCLUDE 'CZT1'
      INCLUDE 'CESTIM'
      INCLUDE 'COUTAU'
      INCLUDE 'CFPLK'
C
C  SAVE INCIDENT SPECIES: IOLD
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,5)
      IOLD=IION
C
      X0=X0+VLXPAR*ZT
      Y0=Y0+VLYPAR*ZT
      Z0=Z0+VLZPAR*ZT
      TIME=TIME+ZT/VELPAR
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=LEARC2(X0,Y0,NRCELL,NPANU,'FOLION 2     ')
      ELSEIF (NLFEM) THEN
        IPOLG=0
      ENDIF
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      MSURF=0
      IF (NLTRA) PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
C
255   CONTINUE
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL UPCUSR(WS,1)
      ENDIF
C
C
C  TEST FOR CORRECT CELL NUMBER AT COLLISION POINT
C  KILL PARTICLE, IF TOO LARGE ROUND OFF ERRORS DURING
C  PARTICLE TRACING
C
      IF (NLTEST) CALL CLLTST(*997)
      IF (ZT.GT.0.D0) THEN
C
C  FLIGHT WITH PARALLEL VELOCITY VELPAR (CM/SEC)
C  PARALLEL DISTANCE ZT (CM)
C  ENERGY RELAXATION CONSTANT TAUE
C
        DUR=ZT/VELPAR
        ENEW=E0*EXP(-DUR/TAUE)+1.5*TI*(1.-EXP(-DUR/TAUE))
        VNEW=RSQDVI(IOLD)*SQRT(ENEW)
C
C  UPDATE ESTIMATORS EIIO,EIPL
        EIIO(NCELL)=EIIO(NCELL)+WEIGHT*(ENEW-E0)
        EIPL(NCELL)=EIPL(NCELL)-WEIGHT*(ENEW-E0)
C
        E0=ENEW
        E0_MEAN=E0
        VEL=VNEW
        RETURN 2
      ENDIF
C
C  POST COLLISION ESTIMATOR
C
C     IF (NCLVI.GT.0) THEN
C       WS=WEIGHT/SIGTOT
C       CALL UPCUSR(WS,2)
C     ENDIF
      RETURN 2
C
997   CALL MASAGE ('ERROR IN FPKCOL,  DETECTED IN SUBR. CLLTST    ')
      CALL MASAGE ('PARTICLE IS KILLED                            ')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,15)
      GOTO 999
C
999   PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      WEIGHT=0.
      CALL LEER(1)
      RETURN 2
      END
C
      SUBROUTINE TORCOL(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  PARTICEL ON TOROIDAL PERIODICITY SURFACE MTSURF, FLIGHT INTO CELL NNTCLL
C  RETURN 1:  NO SURFACE TALLIES, FLIGHT CONTINUES
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'CLGIN'
      INCLUDE 'CLOGAU'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
C
      TIME=TIME+ZT/VEL
C
      IF (NINCZ.EQ.0) THEN
        WRITE (6,*) 'ERROR IN TORCOL, NINCZ ?  '
        WRITE (6,*) 'NINCZ ',NINCZ
        CALL EXIT
      ENDIF
C
      IF (ITYP.LE.2) THEN
        SAVE=VELX
C       SINTOR=SINAL*SIG
        SINTOR=SINAL*FLOAT(NINCZ)
        VELX=COSAL*SAVE+VELZ*SINTOR
        VELZ=-SAVE*SINTOR+COSAL*VELZ
      ENDIF
C  ADVANCE TO NEXT TOROIDAL CELL
      X0=X01-RMTOR
      Y0=Y0+ZT*VELY
      Z01=-Z01
      Z0=Z01
C
      IF (LEVGEO.EQ.3.OR.(LEVGEO.EQ.2.AND.NLPOL)) THEN
        IF (NRCELL.GT.0)
     .  NN=LEARC1(X0,Y0,Z0,IPOLG,NRCELL,NRCELL,.FALSE.,.FALSE.,
     .                           NPANU,'TORCOL 1')
      ENDIF
C
c slmod begin - not tr
c      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,6)
c slmod end
C
      IF (.NOT.NLTOR) THEN
C
        MSURF=0
c slmod begin - tr -new

        NTRSEG=NTRSEG+NINCS
        IF (NTRSEG.LE.0) THEN
          NTRSEG=NTTRA-1
          IF (EIRNSDTOR.GT.1.AND.NACELL.EQ.0) NBLOCK = EIRNSDTOR - 1
        ELSEIF (NTRSEG.GE.NTTRA) THEN
          NTRSEG=1
          IF (EIRNSDTOR.GT.1.AND.NACELL.EQ.0) NBLOCK = 1
        ENDIF

        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,6)

        NLSRFT=.TRUE.

        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(4X,A,I6)') 'TORCOL: NTRSEG=',ntrseg
        ENDIF
c slmod end
        RETURN 1
      ELSE
C
C  NLTOR=TRUE:
C
        MSURF=0
        IF (NRCELL.GT.0) NTCELL=NNTCLL
        ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
        write (6,*) 'torcol ,ntcell,ircell,ipcell,mtsurf,ISTS '
        write (6,*)          ntcell,ircell,ipcell,mtsurf,ists
        IF (ISTS.EQ.0) RETURN 1
        ZT=0.
c slmod begin - tr
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,6)
c slmod end
        RETURN
      ENDIF
C
997   CONTINUE
      END
C
      SUBROUTINE TIMCOL (PR,*,*)
C  COLLISION WITH TIME SURFACE, FIND NEW CO-ORDINATES
C  UPDATE (TIME-) SURFACE TALLIES
C  UPDATE USER SUPPLIED SNAPSHOT ESTIMATED TALLIES (CALL UPNUSR)
C  PUT PARTICLE ONTO CENSUS ARRAYS
C  AND STOP HISTORY
C  RETURN 1: CONTINUE FLIGHT
C  RETURN 2: STOP FLIGHT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'COMNNL'
      INCLUDE 'COMUSR'
      INCLUDE 'CLOGAU'
      INCLUDE 'CLGIN'
      INCLUDE 'CESTIM'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
      DIMENSION RPSTT(NPARTT),IPSTT(MPARTT)
      EQUIVALENCE (RPSTT(1),X0),(IPSTT(1),NPANU)
C
      X0=X0+VELX*TT
      Y0=Y0+VELY*TT
      Z0=Z0+VELZ*TT
      TIME=TIME+TT/VEL
      MSURF=0
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
C
      IPOLG=IPOLGN
      IF (NLTRA) THEN
        PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
        X01=X0+RMTOR
      ENDIF
      X00=X0
      Y00=Y0
      Z00=Z0
      Z01=Z0
      WEIGHT=WEIGHT*PR
C
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,13)
C
C  UPDATE SNAPSHOT ESTIMATORS
      IF (NSNVI.GT.0) CALL UPNUSR
C
c-dpc
      dist=sqrt(x0**2+y0**2)
      if(dist.gt.1e4) then
        write(*,*) 'timcol: ERROR!  dist = ',dist,
     I ISTRA,' (particle more than 100 m from the origin)'
        write(*,*) 'x0,y0,z0,velx,vely,velz ',
     1   x0,y0,z0,velx,vely,velz
        weight=0.
        goto 112
      endif
c-dpc
      IPRNLI=IPRNLI+1
      IPRNLS=IPRNLS+1
      IF (IPRNLS.GE.NPRNLS(ISTRA)) THEN
C  THIS IS THE LAST SCORE FOR THIS STRATUM TO BE STORED
        LGLAST=.TRUE.
      ENDIF
C
C   CENSUS ARRAYS:
C   SAVE LOCATION, WEIGHT AND OTHER PARAMETERS
      DO 100 J=1,NPARTT
        RPART(IPRNLI,J)=RPSTT(J)
100   CONTINUE
      DO 110 J=1,MPARTT
        IPART(IPRNLI,J)=IPSTT(J)
110   CONTINUE
C
112   continue
C  DECIDE: CONTINUE OR STOP TRAJECTORY
      IF ((NTMSTP.GE.0.AND.ITMSTP.GE.NTMSTP).OR.LGLAST) THEN
C
C  DO NOT CONTINUE THIS TRACK
C  UPDATE PARTICLE EFFLUX  ONTO TIME-SURFACE MSURF=NLIM+NSTSI
C  UPDATE ENERGY FLUX ONTO TIME-SURFACE MSURF=NLIM+NSTSI
C
        MSURF=NLIM+NSTSI
        IF (ITYP.EQ.1) THEN
          EOTAT(IATM,MSURF)=EOTAT(IATM,MSURF)+E0*WEIGHT
          POTAT(IATM,MSURF)=POTAT(IATM,MSURF)+WEIGHT
        ELSEIF (ITYP.EQ.2) THEN
          EOTML(IMOL,MSURF)=EOTML(IMOL,MSURF)+E0*WEIGHT
          POTML(IMOL,MSURF)=POTML(IMOL,MSURF)+WEIGHT
        ELSEIF (ITYP.EQ.3) THEN
          EOTIO(IION,MSURF)=EOTIO(IION,MSURF)+E0*WEIGHT
          POTIO(IION,MSURF)=POTIO(IION,MSURF)+WEIGHT
        ENDIF
        ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
        SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        RETURN 2
      ELSE
C  OTHERWISE: RESTORE WEIGHT = WEIGHT/PR, AND CONTINUE
        WEIGHT=WEIGHT/PR
        ITMSTP=ITMSTP+1
        TIME=TIME0
        RETURN 1
      ENDIF
      RETURN
      END
C
      SUBROUTINE ADDCOL (XLI,YLI,ZLI,SG,*,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  SG:    =SIGN OF COSINE OF ANGLE OF INCIDENCE
C  RETURN 1:  NO SURFACE TALLIES, FLIGHT CONTINUES
C  RETURN 2:  LGPART=TRUE:  SURFACE TALLIES,
C                           THEN ABSORBTION, REFLECTION MODEL OR CONTINUATION OF FLIGHT
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMUSR'
      INCLUDE 'CUPD'
      INCLUDE 'CLOGAU'
      INCLUDE 'CADGEO'
      INCLUDE 'CLGIN'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
      INCLUDE 'COUTAU'
      INCLUDE 'CESTIM'
c slmod begin - tr (modified)
      COMMON /TRASHCOM/ nlost
      INTEGER           nlost

      INTEGER CHKVAC
      LOGICAL LTRANS
c slmod end
C   COLLISION WITH ADDITIONAL SURFACE NO. MASURF
C
C  SAVE OLD CO-ORDINATES FOR DIAGNOSTICS
      X0SA=X0
      Y0SA=Y0
      Z0SA=Z0
      MASRFS=MSURF
C  SET NEW CO-ORDINATES ON ADDITIONAL SURFACE. FLIGHT TIME: TL
C  IN CASE OF NLTRA=TRUE, THESE ARE GIVEN IN LOCAL CO-ORDINATE SYSTEM
C  OF TOROIDAL CELL NO ILTOR(MASURF),
C  WHICH MAY OR MAY NOT BE EQUAL NTCELL!
c slmod begin - tr
      IF ( OPTZMOTION.EQ.0.OR.
     .    (OPTZMOTION.EQ.2.AND.NACELL.NE.0).OR.
     .    (OPTZMOTION.EQ.3.AND.NACELL.NE.0.AND.
     .     (X0.GT.62.5.OR.Y0.LT.-61.0)).OR.
     .    (OPTZMOTION.EQ.4.AND.(NACELL.EQ.0.OR.
     .     (X0.LT.62.5.AND.Y0.GT.-61.0)))) Z0=ZLI
      X0=XLI
      Y0=YLI
c
c      X0=XLI
c      Y0=YLI
c      Z0=ZLI
c slmod end
      TIME=TIME+TL/VEL
      MSURF=MASURF
      MRSURF=0
      MPSURF=0
      MTSURF=0
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
      SCOS=SG
      ICOS=SCOS
 
c      if (msurf.EQ.1053) 
c     .  WRITE(0,*) 'MARK: SCOS=',sg

      IF (NLTRA) THEN
C  FIND NEW LOCAL CO-ORDINATES X0,Z0, IN THE CORRECT ZONE
C  FIND TOROIDAL ANGLE PHI
C  AND TEST, IF NTNEW=NTCELL AS IT SHOULD BE
c slmod begin - tr
        IF (ILTOR(MASURF).GT.0.AND.ILTOR(MASURF).NE.NTCELL.AND.
     .      NLTOR) THEN
c
c        IF (ILTOR(MASURF).GT.0.AND.ILTOR(MASURF).NE.NTCELL) THEN
c slmod end
          CALL FZRTOR(X0,Z0,ILTOR(MASURF),XR,PHI,NTNEW,
     .                NLTEST,NTCELL)
          CALL FZRTRI(X0,Z0,NTCELL,XR,PHI,NTCELL)
        ELSE
c slmod begin - tr
          PHI0=PHI
          X02=X01
          Z02=Z01
c slmod end
          PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
        ENDIF
      ENDIF
      IF (NRCELL.GT.0.AND.NRCELL.LE.NR1STM) THEN
        IF (LEVGEO.EQ.3.OR.(LEVGEO.EQ.2.AND.NLPOL))
     .  NN=LEARC1(X0,Y0,Z0,IPOLG,NRCELL,NRCELL,.FALSE.,.FALSE.,
     .                           NPANU,'ADDCOL  ')
      ENDIF
      X00=X0
      Y00=Y0
      Z00=Z0
      Z01=Z0
      IF (NLTRA) X01=X0+RMTOR
C
C
      IWEI=ILSIDE(MASURF)*ICOS
      IF (IWEI.LT.0) GOTO 300
      IF (ILIIN(MASURF).EQ.2) GOTO 400
C
      IF (ILSWCH(MASURF).NE.0) THEN
        NACLLS=NACELL
C  OPERATE A SWITCH
C
C  TURN ON OR OFF THE STANDARD GRID CALCULATION (ONLY ADD. SRF).
        IF (ISWICH(1,MASURF).NE.0) ITIME=ISWICH(1,MASURF)*ICOS
C
C  TURN ON OR OFF MFP AND REACTION RATES: PARTICLE ENTERS VACUUM
        IF (ISWICH(2,MASURF).NE.0) IFPATH=ISWICH(2,MASURF)*ICOS
C
C  TURN ON OR OFF VOLUME AVERAGED TALLIES
        IF (ISWICH(3,MASURF).NE.0) IUPDTE=ISWICH(3,MASURF)*ICOS
C
C  NEW CELL INDEX NACELL, IF PARTICLE IN ADDITIONAL CELL REGION
C  NEW BLOCK INDEX NBLOCK, IF PARTICLE IN STD. MESH REGION
        IF (ISWICH(4,MASURF).NE.0) THEN
          IF (NACELL.GT.0) THEN
c slmod begin - debug - tr - new
            IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10)
     .        WRITE(6,'(3X,A,I4)')
     .          ' ADDCOL: CHANGING ADDITIONAL CELL INDEX ',
     .          NACELL+ISWICH(4,MASURF)*ICOS*ILACLL(MASURF)

            IF     (OPTVAC.GT.0.AND.ILACLL(MASURF).EQ.999) THEN
              IF (NLTRA) THEN
                IF (NLTOR) THEN
                ELSE
                  PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
                  TOR=(ATAN2(Z0,X0)+DBLE(NTRSEG-1)*PHISEG)*RADDEG
c...              Special case for (almost) full toroidal grid:
                  IF (NTRSEG.EQ.1.AND.TOR.LT.0.0D0) TOR = TOR + 360.0D0
                ENDIF
              ELSE
                TOR=Z0
              ENDIF
              NACELL=CHKVAC(-1,NACELL,X0,Y0,TOR,3,NPANU)
c              NACELL = CHKVAC(-1,NACELL,X0,Y0,Z0,3,NPANU)
              IF     (NACELL.EQ.-1) THEN
                WRITE(6,*) 'ADDITIONAL CELL NOT FOUND'
                IWEI=-10
                NACELL=-2
                GOTO 300
              ELSEIF (NACELL.EQ.-2) THEN
                WRITE(6,*) 'MORE THAN ONE ADDITIONAL CELL FOUND'
                IWEI=-10
                NACELL=-2
                GOTO 300
              ENDIF
            ELSEIF (OPTVAC.GT.0.AND.ILACLL(MASURF).EQ.997) THEN
              IF (NLTRA) THEN
                IF (NLTOR) THEN
                ELSE
                  PHISEG=2.0D0*PIA/DBLE(NTTRA-1)
                  TOR=(ATAN2(Z0,X0)+DBLE(NTRSEG-1)*PHISEG)*RADDEG
c...              Special case for (almost) full toroidal grid:
                  IF (NTRSEG.EQ.1.AND.TOR.LT.0.0D0) TOR = TOR + 360.0D0
                ENDIF
              ELSE
                TOR=Z0
              ENDIF
              OLDNACELL=NACELL
              NACELL=CHKVAC(-1,NACELL,X0,Y0,TOR,3,NPANU)
c              NACELL=CHKVAC(-1,NACELL,X0,Y0,Z0,3,NPANU)
              IF     (NACELL.EQ.-1) THEN
                IF (OLDNACELL.EQ.1) THEN
                  WRITE(6,*) 'UNABLE TO FIND ADDITIONAL CELL'
c                  optvac = 2 
c                  NACELL = CHKVAC(-1,NACELL,X0,Y0,Z0)
c                  optvac = 1
                  IWEI=-10
                  NACELL=-2
                  GOTO 300
                ELSE
c ASSUME THAT THE PARICLE IS ENTERING ADDITIONAL CELL 1
                  NACELL=1
                ENDIF
              ELSEIF (NACELL.EQ.-2) THEN
                WRITE(6,*) 'MORE THAN ONE ADDITIONAL CELL IDENTIFIED'
                IWEI=-10
                NACELL=-2
                GOTO 300
              ENDIF
            ELSEIF (ILACLL(MASURF).EQ.998) THEN
c...          Transparent toroidal end surface:
              IF     (ASC3DMODE.EQ.1) THEN
                STOP 'NOT DONE - A'
              ELSEIF (ASC3DMODE.EQ.2) THEN
c...            Less general 3D vacuum mesh:
                IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,98,6)
                IF     (MASURF.EQ.NLIMI-1) THEN
                  IF (NLTRA) THEN
                    IF (NLTOR) THEN
                    ELSE
                      Z0=X0*P1(3,NLIMI)/P1(1,NLIMI)
                      PHI=ATAN2(Z0,X0)
                      NTRSEG=ILTOR(NLIMI)
                    ENDIF
                  ELSE
                    Z0=ZAA
                  ENDIF
                  MASURF=NLIMI
c...              THIS WILL FAIL FOR TOROIDALLY CONTINUOUS PRESSURE GAUGES! FIX!
                  IF (NACELL.GT.1) NACELL=NACELL+(ASCNCUT-1)*ASCNCELL
                ELSEIF (MASURF.EQ.NLIMI) THEN
                  IF (NLTRA) THEN
                    IF (NLTOR) THEN
                    ELSE
                      Z0=X0*P1(3,NLIMI-1)/P1(1,NLIMI-1)
                      PHI=ATAN2(Z0,X0)
                      NTRSEG=ILTOR(NLIMI-1)
                    ENDIF
                  ELSE
                    Z0=0.0D0
                  ENDIF
                  MASURF=NLIMI-1
                  IF (NACELL.GT.1) NACELL=NACELL-(ASCNCUT-1)*ASCNCELL
                ELSE
                  WRITE(0,*) 'CRAP - ON ADDITIONAL MESH - 2D=2'
                  WRITE(0,*) 'Z0,ZAA=',z0,zaa
                  WRITE(0,*) 'MASURF,EBGKI=',masurf,ebgki,npanu
                  STOP 
                ENDIF

c                IF (NLTRA) THEN
c                  IF (NLTOR) THEN
c                  ELSE
c                    WRITE(0,*) 'NEED DEV'
c                  ENDIF
c                ELSE
c...              Cylindrical approximation:
c                  IF     (Z0.EQ.ZAA  .AND.MASURF.EQ.NLIMI) THEN
c                    Z0=0.0D0
c..THIS MEANS THE BGK GRID MUST EXIST?  ASSUMES THAT THE TOROIDAL BOUNDARYS ARE THE NEXT 2 SURFACES!  NOT GARUNTEED!
c                    MASURF=NLIMI-1
c... THIS WILL FAIL FOR TOROIDALLY CONTINUOUS PRESSURE GAUGES! FIX!
c                    IF (NACELL.GT.1) NACELL=NACELL-(ASCNCUT-1)*ASCNCELL
c                  ELSEIF (Z0.EQ.0.0D0.AND.MASURF.EQ.NLIMI-1) THEN
c                    Z0=ZAA
c                    MASURF=NLIMI
c... THIS WILL FAIL FOR TOROIDALLY CONTINUOUS PRESSURE GAUGES!
c                    IF (NACELL.GT.1) NACELL=NACELL+(ASCNCUT-1)*ASCNCELL
c                  ELSE
c                    WRITE(0,*) 'CRAP - ON ADDITIONAL MESH - 2D=2'
c                    WRITE(0,*) 'Z0,ZAA=',z0,zaa
c                    WRITE(0,*) 'MASURF,EBGKI=',masurf,ebgki,npanu
c                    STOP 
c                  ENDIF
c                ENDIF

              ELSE
                IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,98,6)
                IF (NLTRA) THEN
                  IF (NLTOR) THEN
                  ELSE
                    IF     (MASURF.EQ.NLIMI-1) THEN
                      Z0=X0*P1(3,NLIMI)/P1(1,NLIMI)
                      PHI=ATAN2(Z0,X0)
                      NTRSEG=ILTOR(NLIMI)
                      MASURF=NLIMI
                    ELSEIF (MASURF.EQ.NLIMI) THEN
                      Z0=X0*P1(3,NLIMI-1)/P1(1,NLIMI-1)
                      PHI=ATAN2(Z0,X0)
                      NTRSEG=ILTOR(NLIMI-1)
                      MASURF=NLIMI-1
                    ELSE
                      WRITE(0,*) 'CRAP - ON ADDITIONAL MESH - zuper!'
                      WRITE(0,*) 'MASURF,EBGKI=',masurf,ebgki
                    ENDIF
                  ENDIF
                ELSE
                  IF     (Z0.EQ.ZAA  .AND.MASURF.EQ.NLIMI) THEN
                    Z0=0.0D0
                    MASURF=NLIMI-1
                  ELSEIF (Z0.EQ.0.0D0.AND.MASURF.EQ.NLIMI-1) THEN
                    Z0=ZAA
                    MASURF=NLIMI
                  ELSE
                    WRITE(0,*) 'CRAP - ON ADDITIONAL MESH'
                    WRITE(0,*) 'Z0,ZAA=',z0,zaa
                    WRITE(0,*) 'MASURF,EBGKI=',masurf,ebgki
                    STOP 
                  ENDIF
                ENDIF
              ENDIF
              MSURF=MASURF
              X00=X0
              Y00=Y0
              Z00=Z0
              Z01=Z0
              IF (NLTRA) X01=X0+RMTOR
            ELSE
              NACELL=NACELL+ISWICH(4,MASURF)*ICOS*ILACLL(MASURF)
            ENDIF
c
c            NACELL=NACELL+ISWICH(4,MASURF)*ICOS*ILACLL(MASURF)
c slmod end
            IF (NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
              IWEI=-10
              GOTO 300
            ENDIF
          ELSEIF (NACELL.EQ.0) THEN
c slmod begin - debug - tr
            IF (ILBLCK(MASURF).EQ.998) THEN
c...          Transparent toroidal end surface:

              IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,98,6)

              IF (NLTRA) THEN  
c...            Toroidal approximation:
                IF (NLTOR) THEN
                ELSE
                  IF     (MASURF.EQ.NLIMI) THEN
                    Z0=X0*P1(3,NLIMI-1)/P1(1,NLIMI-1)
                    PHI=ATAN2(Z0,X0)
                    NTRSEG=ILTOR(NLIMI-1)
                    MASURF=NLIMI-1
                    IF (EIRNSDTOR.GT.1) NBLOCK=1
                  ELSEIF (MASURF.EQ.NLIMI-1) THEN
                    Z0=X0*P1(3,NLIMI)/P1(1,NLIMI)
c                      PHI=MOD(PHI0-ATAN2(Z02,X02)+
c     .                             ATAN2(Z0,(RMTOR+X0)),PI2A)
                    PHI=ATAN2(Z0,X0)
                    NTRSEG=ILTOR(NLIMI)
                    MASURF=NLIMI
                    IF (EIRNSDTOR.GT.1) NBLOCK=EIRNSDTOR-1
                  ELSE
                    WRITE(0,*) 'CRAP - ON STANDARD MESH'
                    WRITE(0,*) 'Z0,ZAA=',z0,zaa
                    WRITE(0,*) 'MASURF,EBGKI=',masurf,ebgki
                    STOP 
                  ENDIF
                ENDIF

	      ELSE   
c...            Cylindrical approximation:
                IF     (DABS(Z0-ZAA).LT.EPS10.AND.MASURF.EQ.NLIMI) THEN
                  Z0=0.0D0
c..THIS MEANS THE BGK GRID MUST EXIST?
                  MASURF=NLIMI-1
                  IF (EIRNSDTOR.GT.1) THEN
                    NBLOCK=1
                  ENDIF
                ELSEIF (DABS(Z0).LT.EPS10.AND.MASURF.EQ.NLIMI-1) THEN
                  Z0=ZAA
                  MASURF=NLIMI
                  IF (EIRNSDTOR.GT.1) THEN
                    NBLOCK=EIRNSDTOR-1
                  ENDIF
                ELSE
                  WRITE(0,*) 'CRAP - ON STANDARD MESH (CYL)'
                  WRITE(0,*) 'Z0,ZAA=',z0,zaa
                  WRITE(0,*) 'MASURF,EBGKI=',masurf,ebgki
                  STOP 
                ENDIF
              ENDIF

              IF (NBLOCK.GT.NBMLT.OR.NBLOCK.LT.1) THEN
                IWEI=-10
                GOTO 300
              ENDIF

              MSURF=MASURF
              X00=X0
              Y00=Y0
              Z00=Z0
              Z01=Z0
              IF (NLTRA) X01=X0+RMTOR
            ELSE
              NBLOCK=NBLOCK+ILBLCK(MASURF)*ICOS*ISWICH(4,MASURF)
              IF (.NOT.NLMLT.OR.(NBLOCK.GT.NBMLT.OR.NBLOCK.LT.1)) THEN
                WRITE(6,*) 'MARK: NBLOCK CRAP OUT'
                IWEI=-10
                GOTO 300
              ENDIF
            ENDIF
c
c            NBLOCK=NBLOCK+ILBLCK(MASURF)*ICOS*ISWICH(4,MASURF)
c            IF (.NOT.NLMLT.OR.(NBLOCK.GT.NBMLT.OR.NBLOCK.LT.1)) THEN
c               IWEI=-10
c              GOTO 300
c            ENDIF
c slmod end
          ENDIF
          NBLCKA=NSTRD*(NBLOCK-1)+NACELL
C
C  ENTRANCE INTO STANDARD MESH, INTO BLOCK NBLOCK=ILBLCK
C  OR
C  EXIT FROM STANDARD MESH, INTO CELL NACELL=ILACLL
        ELSEIF (ISWICH(5,MASURF).NE.0) THEN
          IF (NACELL.EQ.0) THEN
C  SET CELL INDEX EQUAL TO ILACLL
            NACELL=ILACLL(MASURF)
            NBLOCK=NBMLTP
C
            NRCELL=0
            NPCELL=1
            NTCELL=1
            IF (.NOT.NLADD.OR.NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
              IWEI=-10
              GOTO 300
            ENDIF
          ELSEIF (NACELL.GT.0) THEN
            STOP 'MARK: STOP: NOT ALLOWED TO ENTER STANDARD MESH'
C  FIND  NRCELL,IPOLG,NPCELL,NTCELL IN STANDARD MESH, BLOCK ILBLCK
            NACELL=0
            NBLOCK=ILBLCK(MASURF)
C
            IAN=1
            IEN=NR1STM
            NRCELL=LEARC1(X0,Y0,Z0,IPOLG,1,NR1STM,.FALSE.,.FALSE.,NPANU,
     .                   'ADDCOL      ')
            IF (NLTOR) THEN
              IF (NLTRZ) THEN
                NTCELL=LEARCA(Z0,ZSURF,1,NT3RD,1,'ADDCOL   ')
              ELSEIF (NLTRA) THEN
                NTCELL=LEARCA(PHI,ZSURF,1,NT3RD,1,'ADDCOL   ')
              ENDIF
            ELSE
              NTCELL=1
            ENDIF
            IF (NLPOL) THEN
              IF (LEVGEO.EQ.1) THEN
                NPCELL=LEARCA(Y0,PSURF,1,NP2ND,1,'ADDCOL')
              ELSEIF (LEVGEO.EQ.2) THEN
                IF (NLCRC) THEN
                  WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
                  NPCELL=LEARCA(WINK,PSURF,1,NP2ND,1,'ADDCOL')
                ELSE
                  NPCELL=LEARC2(X0,Y0,NRCELL,NPANU,'ADDCOL')
                ENDIF
              ELSEIF (LEVGEO.EQ.3) THEN
                NPCELL=IPOLG
              ELSE
                WRITE (6,*) 'ERROR EXIT FROM ADDCOL. NLPOL ',LEVGEO
                CALL EXIT
              ENDIF
            ELSE
              NPCELL=1
            ENDIF
          ENDIF
          NBLCKA=NSTRD*(NBLOCK-1)+NACELL
C
        ELSEIF (ISWICH(6,MASURF).NE.0) THEN
          IF (NACELL.EQ.0) THEN
C  SET CELL INDEX EQUAL TO NBLOCK+ILACLL
            NACELL=NBLOCK+ICOS*ISWICH(6,MASURF)*ILACLL(MASURF)
            NBLOCK=NBMLTP
C
            NRCELL=0
            NPCELL=1
            NTCELL=1
            IF (.NOT.NLADD.OR.NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
              IWEI=-10
              GOTO 300
            ENDIF
          ELSEIF (NACELL.GT.0) THEN
C  FIND  NRCELL,IPOLG,NPCELL,NTCELL IN STANDARD MESH, BLOCK ILBLCK
            NBLOCK=NACELL+ICOS*ISWICH(6,MASURF)*ILBLCK(MASURF)
            NACELL=0
C
            IAN=1
            IEN=NR1STM
            NRCELL=LEARC1(X0,Y0,Z0,IPOLG,1,NR1STM,.FALSE.,.FALSE.,NPANU,
     .                   'ADDCOL7     ')
            IF (NLTOR) THEN
              IF (NLTRZ) THEN
                NTCELL=LEARCA(Z0,ZSURF,1,NT3RD,1,'ADDCOL    ')
              ELSEIF (NLTRA) THEN
                NTCELL=LEARCA(PHI,ZSURF,1,NT3RD,1,'ADDCOL   ')
              ENDIF
            ELSE
              NTCELL=1
            ENDIF
            IF (NLPOL) THEN
              IF (LEVGEO.EQ.1) THEN
                NPCELL=LEARCA(Y0,PSURF,1,NP2ND,1,'ADDCOL')
              ELSEIF (LEVGEO.EQ.2) THEN
                IF (NLCRC) THEN
                  WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
                  NPCELL=LEARCA(WINK,PSURF,1,NP2ND,1,'ADDCOL')
                ELSE
                  NPCELL=LEARC2(X0,Y0,NRCELL,NPANU,'ADDCOL')
                ENDIF
              ELSEIF (LEVGEO.EQ.3) THEN
                NPCELL=IPOLG
              ELSE
                WRITE (6,*) 'ERROR EXIT FROM ADDCOL. NLPOL ',LEVGEO
                CALL EXIT
              ENDIF
            ELSE
              NPCELL=1
            ENDIF
          ENDIF
          NBLCKA=NSTRD*(NBLOCK-1)+NACELL
C
        ENDIF
C
      ENDIF
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,6)
C
c slmod begin - not tr
c      IF (npanu.EQ.425) THEN
c        WRITE(80,*) 'MASURF,MSURF=',MASURF,MSURF,ILIIN(MASURF)
c      ENDIF
c slmod end
      IF (ILIIN(MASURF).LT.0) THEN
        IF (ILIIN(MASURF).EQ.-1) RETURN 1
        RETURN 2
      ENDIF
C
C  ILIIN(MASURF) .EQ. 1, OR ILIIN(MASURF) .EQ.3
C  PREPARE REFLECTION, I.E. SET OUTER NORMAL
C
      IF (JUMLIM(MASURF).GT.0) GOTO 200
      GOTO 100
C
      ENTRY ADDNOR (X0E,Y0E,Z0E,SCOSE,MSURFE,*,*)
      X0=X0E
      Y0=Y0E
      Z0=Z0E
      SCOS=SCOSE
      MSURF=MSURFE
      IF (JUMLIM(MASURF).GT.0) GOTO 200
C
100   CONTINUE
      CRTX=A1LM(MASURF)+ALM(MASURF)*X0+A7LM(MASURF)*Y0+A8LM(MASURF)*Z0
      CRTY=A2LM(MASURF)+A7LM(MASURF)*X0+BLM(MASURF)*Y0+A9LM(MASURF)*Z0
      CRTZ=A3LM(MASURF)+A8LM(MASURF)*X0+A9LM(MASURF)*Y0+CLM(MASURF)*Z0
      CNORM=1./SQRT(CRTX*CRTX+CRTY*CRTY+CRTZ*CRTZ)*SCOS
      CRTX=CRTX*CNORM
      CRTY=CRTY*CNORM
      CRTZ=CRTZ*CNORM
      RETURN 2
200   CONTINUE
      CRTX=A1LM(MASURF)*SCOS
      CRTY=A2LM(MASURF)*SCOS
      CRTZ=A3LM(MASURF)*SCOS
      RETURN 2
C
C
300   CONTINUE
      IF (IWEI.EQ.-1) THEN
C  PARTICLE HAS HIT A SURFACE FROM AN ABSORBING SIDE
C  UPDATE FLUXES (DO NOT SET WEIGHT=0.D0) AND ABSORB PARTICLE
        IF (NLTRC) THEN
          CALL CHCTRC(X0,Y0,Z0,16,6)
          WRITE (6,*) 'ABSORB PARTICLE: NPANU ',NPANU
        ENDIF
        SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        LGPART=.FALSE.
        RETURN 2
      ELSEIF (IWEI.EQ.-2) THEN
        PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
        ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
C  KILL THIS PARTICLE BECAUSE OF TOO LARGE ROUND-OFF ERRORS DURING
C  THE PARTICLE-TRACING. DO NOT UPDATE FLUXES (SET WEIGHT=0.D0)
        IF (NLTRC) THEN
          CALL CHCTRC(X0,Y0,Z0,16,15)
        ENDIF
        WRITE (6,*) 'ERROR DETECTED IN SUBR. ADDCOL '
        WRITE (6,*) 'PARTICLE COMES FROM WRONG SIDE '
        CALL MASJ1 ('NPANU=  ',NPANU)
        CALL MASJ1 ('MASRF NW',MASURF)
        CALL MASJ1 ('MASRF OD',MASRFS)
        CALL MASR3 ('X0,Y0,Z0 (NEW)          ',X0,Y0,Z0)
        CALL MASR3 ('X0,Y0,Z0 (OLD)          ',X0SA,Y0SA,Z0SA)
        CALL MASR3 ('VELX,VELY,VELZ          ',VELX,VELY,VELZ)
        CALL MASR2 ('WEIGHT,E0       ',WEIGHT,E0)
        SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        WEIGHT=0.
        LGPART=.FALSE.
c slmod begin - tr - new
        NLOST=NLOST+1
c slmod end
        RETURN 2
      ELSEIF (IWEI.EQ.-3) THEN
C  SURFACE IS NOT SEEN BY THE PARTICLE BECAUSE OF ILSIDE OPTION
C  I.E. SURFACE IS TRANSPARENT FROM THIS SIDE
C  ACTS AS ILIIN=0 OPTION (NO SURFACE TALLIES, NO SWITCHES)
        RETURN 1
      ELSEIF (IWEI.EQ.-10) THEN
C  KILL THIS PARTICLE BECAUSE CELL NUMBER OF RANGE DUE TO SWITCHING
C  DO NOT UPDATE FLUXES (SET WEIGHT=0.D0)
        PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
        ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
        IF (NLTRC) THEN
          CALL CHCTRC(X0,Y0,Z0,16,15)
        ENDIF
        WRITE (6,*) 'ERROR DETECTED IN SUBR. ADDCOL '
        WRITE (6,*) 'SOME CELL INDEX OUT OF RANGE '
        CALL MASJ1 ('NPANU=  ',NPANU)
        WRITE (6,*) 'NLMLT,NLADD ',NLMLT,NLADD
        WRITE (6,*) 'NBMLT,NRADD ',NBMLT,NRADD
        CALL MASJ1 ('MASRF NW',MASURF)
        CALL MASJ1 ('MASRF OD',MASRFS)
        CALL MASJ1 ('NACL NEW',NACELL)
        CALL MASJ1 ('NACL OLD',NACLLS)
        CALL MASJ1 ('NBLOCK  ',NBLOCK)
        CALL MASR3 ('X0,Y0,Z0 (NEW)          ',X0,Y0,Z0)
        CALL MASR3 ('X0,Y0,Z0 (OLD)          ',X0SA,Y0SA,Z0SA)
        CALL MASR3 ('VELX,VELY,VELZ          ',VELX,VELY,VELZ)
        CALL MASR2 ('WEIGHT,E0       ',WEIGHT,E0)
        SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        WEIGHT=0.
        LGPART=.FALSE.
c slmod begin - tr
        NLOST=NLOST+1
c slmod end
        RETURN 2
      ENDIF
C
400   CONTINUE
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,6)
      RETURN 2
      END
C
      SUBROUTINE ESCAPE(PR,SG,*,*,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMUSR'
      INCLUDE 'CSPEZ'
      INCLUDE 'CESTIM'
      INCLUDE 'COMSPL'
      INCLUDE 'CLGIN'
      INCLUDE 'CRAND'
      INCLUDE 'CCONA'
      INCLUDE 'CZT1'
      INCLUDE 'COUTAU'
      INCLUDE 'CGRID'
      INCLUDE 'CLOGAU'
      LOGICAL NLSPUT,LTRANS
      DIMENSION DIWL(NPLS),VPWL(NPLS)
c slmod begin - debug - tr
      COMMON /DELTAE0COM/ DELTAE0,NDELTAE0
      REAL*8              DELTAE0,NDELTAE0

      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: NO DEBUG INFROMATION IN '//
     .             'THIS VERSION'


c slmod end
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
        IF (.NOT.LGPART) THEN
          WRITE (6,*) 'ERROR AT PERIODICITY SURFACE, LGPART=FALSE '
          RETURN
        ENDIF
        IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,0,9)
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 1'
c slmod end
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
C  FOR TRANSPARENT SURFACES: ONE SIDED FLUX, POSITIVE COMPONENT
C
c slmod begin - debug - not tr
c        IF (DEBUGOPT.NE.0) THEN
c          WRITE(0,*) 'ESCAPE: BRANCH =',
c     .      (ILIIN(MSURF).LT.0).AND.(SG.LT.0.D0),MSURF,SG
c          WRITE(6,*) 'ESCAPE: BRANCH =',
c     .      (ILIIN(MSURF).LT.0).AND.(SG.LT.0.D0),MSURF,SG
c        ENDIF
c      WRITE(0,*) 'POTAT ?:',msurf,sg
c...bug
c...I can't recall why I wanted to do this.  I think it was to 
c   calculate the net flux of neutrals through a transparent surface,
c   and limit selection based on orientation to surfaces with
c   ILIIN.EQ.-1.  Not sure however.
c      IF ((ILIIN(MSURF).GT.-2).AND.(SG.LT.0.D0)) GOTO 10
c
c      IF (msurf.EQ.1119) THEN
c        WRITE(0,*) '1119:',ILIIN(MSURF),sg,npanu
c      ENDIF

      IF ((ILIIN(MSURF).LT.0).AND.(SG.LT.0.D0)) GOTO 10
c slmod end
C
C  UPDATE PARTICLE EFFLUX  ONTO SURFACE MSURF
C  UPDATE ENERGY FLUX ONTO SURFACE MSURF
C
C  SPATIAL RESOLUTION ON NON DEFAULT STANDARD SURFACE?
c slmod begin - juelich - not tr (already fixed)
      IF (MSURF.GT.NLIM.AND.LEVGEO.LT.4.AND.NLMPGS.GT.NLIMPS) THEN
c
c      IF (MSURF.GT.NLIM.AND.LEVGEO.LT.4) THEN
c slmod end
        ISTS=MSURF-NLIM
        IF (INUMP(ISTS,1).NE.0) MSURFG=NPCELL+(NTCELL-1)*NP2T3
        IF (INUMP(ISTS,2).NE.0) MSURFG=NRCELL+(NTCELL-1)*NR1P2
        IF (INUMP(ISTS,3).NE.0) MSURFG=NRCELL+(NPCELL-1)*NR1P2
        MSURFG=NLIM+NSTSI+MSURFG+(ISTS-1)*NGITT
        FLX=FLXOUT(MSURFG)
c slmod begin - debug - not tr
c        IF (DEBUGOPT.NE.0) THEN
c          WRITE(0,*) 'ESCAPE: ASSIGNING MSURFG 01 =',msurfg
c          WRITE(6,*) 'ESCAPE: ASSIGNING MSURFG 01 =',msurfg
c        ENDIF
c
c        WRITE(6,*) 'MARK: SPATIAL RESOLUTION - MSURFG'
c        WRITE(6,'(8I6)') INUMP(ISTS,1),ISTS,MSURFG,NPCELL,NTCELL,NP2T3,
c     .                   NGITT,MSURF
c        WRITE(6,'(8I6)') INUMP(ISTS,2),ISTS,MSURFG,NRCELL,NTCELL,NR1P2,
c     .                   NGITT,MSURF
c        WRITE(6,'(8I6)') INUMP(ISTS,3),ISTS,MSURFG,NRCELL,NPCELL,NR1P2,
c     .                   NGITT,MSURF
c slmod end
      ELSE
        MSURFG=0
        FLX=FLXOUT(MSURF)
c slmod begin - debug - not tr
c        IF (DEBUGOPT.NE.0) THEN
c          WRITE(0,*) 'ESCAPE: ASSIGNING MSURFG 01 =',msurfg
c          WRITE(6,*) 'ESCAPE: ASSIGNING MSURFG 01 =',msurfg
c        ENDIF
c slmod end
      ENDIF
C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: B E0= ',E0,ityp,msurfg
c slmod end
      IF (ITYP.EQ.1) THEN
        EOTAT(IATM,MSURF)=EOTAT(IATM,MSURF)+E0*WPR
        POTAT(IATM,MSURF)=POTAT(IATM,MSURF)+WPR
        FMASS=DBLE(NMASSA(IATM))
        FCHAR=DBLE(NCHARA(IATM))
      ELSEIF (ITYP.EQ.2) THEN
        EOTML(IMOL,MSURF)=EOTML(IMOL,MSURF)+E0*WPR
        POTML(IMOL,MSURF)=POTML(IMOL,MSURF)+WPR
c        IF (msurf.EQ.1053) 
c     .    WRITE(0,'(A,I6,3F10.5,I6,1P,E12.4)') 
c     .      'MARK: WPR:',npanu,pr,sg,wpr,imol,potml(imol,msurf)
        FMASS=DBLE(NMASSM(IMOL))
        FCHAR=DBLE(NCHARM(IMOL))
      ELSEIF (ITYP.EQ.3) THEN
        IF (ILIIN(MSURF).GT.0) THEN
C  ACCOUNT FOR ELECTROSTATIC SHEATH AT SURFACE FOR TEST IONS
          IF (FSHEAT(MSURF).LE.0.D0) THEN
            GAMMA=0.
            CUR=0.
            IC=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
            IF (.NOT.LGVAC(IC,NPLSI+1)) THEN
              TEWL=TEIN(IC)
              DO 30 IP=1,NPLSI
                VPWL(IP)=SQRT(VXIN(IP,IC)**2+VYIN(IP,IC)**2+
     .                        VZIN(IP,IC)**2)
                DIWL(IP)=DIIN(IP,IC)
30            CONTINUE
              ESHET=NCHRGI(IION)*SHEATH(TEWL,DIWL,VPWL,
     .                                  NCHRGP,GAMMA,CUR,NPLSI,MSURF)
            ENDIF
          ELSE
            ESHET=FSHEAT(MSURF)*TEWL*NCHRGI(IION)
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
        EOTIO(IION,MSURF)=EOTIO(IION,MSURF)+E0*WPR
        POTIO(IION,MSURF)=POTIO(IION,MSURF)+WPR
        FMASS=DBLE(NMASSI(IION))
        FCHAR=DBLE(NCHARI(IION))
      ENDIF
      IF (MSURFG.GT.0) THEN
        IF (ITYP.EQ.1) THEN
          EOTAT(IATM,MSURFG)=EOTAT(IATM,MSURFG)+E0*WPR
          POTAT(IATM,MSURFG)=POTAT(IATM,MSURFG)+WPR
        ELSEIF (ITYP.EQ.2) THEN
          EOTML(IMOL,MSURFG)=EOTML(IMOL,MSURFG)+E0*WPR
          POTML(IMOL,MSURFG)=POTML(IMOL,MSURFG)+WPR
        ELSEIF (ITYP.EQ.3) THEN
          EOTIO(IION,MSURFG)=EOTIO(IION,MSURFG)+E0*WPR
          POTIO(IION,MSURFG)=POTIO(IION,MSURFG)+WPR
        ENDIF
      ENDIF
      ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
C
10    CONTINUE
C
C  ADDITIONAL INCIDENT SURFACE TALLIES
      IF (NADSI.GE.1) CALL UPSUSR (WPR,1)
c slmod begin - debug - tr
c        IF (DEBUGOPT.NE.0) THEN
c          WRITE(0,*) 'ESCAPE:           MSURFG 02 =',msurfg
c          WRITE(6,*) 'ESCAPE:           MSURFG 02 =',msurfg
c        ENDIF
c slmod end
C
C  STOP TRAJECTORY, FOR SOME REASON IN SUBROUTINE ADDCOL OR STDCOL
C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN ?',e0
c slmod end
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
      WGHTSS=0.
      YIELD1=0.
      YIELD2=0.
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
        WGHTSS=WPR*YIELD2
C
C  UPDATE SPUTTERED FLUX IF AVAILABLE. SORTED BY INCIDENT PARTICLE TYPE
C
        IF (NLSPUT) THEN
          IF (ITYP.EQ.1) THEN
            SPTAT(IATM,MSURF)=SPTAT(IATM,MSURF)+WGHTSP+WGHTSS
          ELSEIF (ITYP.EQ.2) THEN
            SPTML(IMOL,MSURF)=SPTML(IMOL,MSURF)+WGHTSP+WGHTSS
          ELSEIF (ITYP.EQ.3) THEN
            SPTIO(IION,MSURF)=SPTIO(IION,MSURF)+WGHTSP+WGHTSS
C         ELSEIF (ITYP.EQ.4) THEN
C           SPTPL(IPLS,MSURF)=SPTPL(IPLS,MSURF)+WGHTSP+WGHTSS
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (ITYP.EQ.1) THEN
              SPTAT(IATM,MSURFG)=SPTAT(IATM,MSURFG)+WGHTSP+WGHTSS
            ELSEIF (ITYP.EQ.2) THEN
              SPTML(IMOL,MSURFG)=SPTML(IMOL,MSURFG)+WGHTSP+WGHTSS
            ELSEIF (ITYP.EQ.3) THEN
              SPTIO(IION,MSURFG)=SPTIO(IION,MSURFG)+WGHTSP+WGHTSS
C           ELSEIF (ITYP.EQ.4) THEN
C             SPTPL(IPLS,MSURFG)=SPTPL(IPLS,MSURFG)+WGHTSP+WGHTSS
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
c slmod begin - debug - not tr
c        IF (DEBUGOPT.NE.0) THEN
c          WRITE(0,*) 'ESCAPE:           MSURFG 03 =',msurfg
c          WRITE(6,*) 'ESCAPE:           MSURFG 03 =',msurfg
c        ENDIF
c slmod end
      LTRANS=.FALSE.
      IF (TRANSP(1,MSURF).GT.0.D0.OR.TRANSP(2,MSURF).GT.0.D0) THEN
C
C  AT THIS POINT: ILIIN(MSURF).GT.0
C
        IF (SG.GT.0) ISG=1
        IF (SG.LT.0) ISG=2
        LTRANS=RANF_EIRENE( ).LE.TRANSP(ISG,MSURF)
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
            IF (IDIM.EQ.2) NPCELL=NPCELL+NINCY
            IF (IDIM.EQ.3) NTCELL=NTCELL+NINCZ
          ENDIF
          IF (NLTRC) THEN
            CALL LEER(1)
            WRITE (6,*) 'SURFACE MSURF= ',MS,' IS MADE TRANSPARENT'
          ENDIF
        ENDIF
C
      ENDIF
c slmod begin - not tr
c      WRITE(6,*) 'MARK: NPANU,NA,NR,NP,LTRANS=',
c     .                  NPANU,NACELL,NRCELL,NPCELL,LTRANS
c slmod end
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
      IF (ILIIN(MSURF).EQ.2) THEN
        SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WPR
        WEIGHT=0.D0
        LGPART=.FALSE.
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN - ABSORBING'
c slmod end
        RETURN
      ENDIF
C
C   ..............................................
C   .                                            .
C   .   MIRROR, OR SEMI-TRANSPARENT SURFACE      .
C   .   REEMITTED FLUX=INCOMING FLUX AND RETURN  .
C   ..............................................
C
c slmod begin - debug - not tr
c        IF (DEBUGOPT.NE.0) THEN
c          WRITE(0,*) 'ESCAPE:           MSURFG 04 =',msurfg
c          WRITE(6,*) 'ESCAPE:           MSURFG 04 =',msurfg
c        ENDIF
c slmod end
      IF (LTRANS.OR.ILIIN(MSURF).EQ.3) THEN
C
C ITOLD=ITNEW=ITYP
C
        IF (ITYP.EQ.1) THEN
          ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WPR
          PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WPR
        ELSEIF (ITYP.EQ.2) THEN
          ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WPR
          PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WPR
        ELSEIF (ITYP.EQ.3) THEN
          ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WPR
          PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WPR
        ENDIF
        IF (MSURFG.GT.0) THEN
          IF (ITYP.EQ.1) THEN
            ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WPR
            PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WPR
          ELSEIF (ITYP.EQ.2) THEN
            ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WPR
            PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WPR
          ELSEIF (ITYP.EQ.3) THEN
            ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WPR
            PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WPR
          ENDIF
        ENDIF
C
C  EITHER: SEMI-TRANSPARENT SURFACE
C
        IF (LTRANS) THEN
          IF (NADSI.GE.1) CALL UPSUSR (WPR,2)
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 2'
c slmod end
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
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 3'
c slmod end
          RETURN 1
        ENDIF
      ENDIF
C
C   .........................
C   .                       .
C   .  TRANSPARENT SURFACE  .
C   .........................
C
c slmod begin - debug - not tr
c        IF (DEBUGOPT.NE.0) THEN
c          WRITE(0,*) 'ESCAPE:           MSURFG 05 =',msurfg
c          WRITE(6,*) 'ESCAPE:           MSURFG 05 =',msurfg
c        ENDIF
c slmod end
      IF (ILIIN(MSURF).LT.0) THEN
C
C  ONE SIDED FLUX: NEGATIVE COMPONENT
C
        IF (SG.GT.0.D0) GOTO 90
C
C ITOLD=ITNEW=ITYP
C
        IF (ITYP.EQ.1) THEN
          ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WPR
          PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WPR
        ELSEIF (ITYP.EQ.2) THEN
          ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WPR
          PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WPR
        ELSEIF (ITYP.EQ.3) THEN
          ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WPR
          PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WPR
        ENDIF
        IF (MSURFG.GT.0) THEN
          IF (ITYP.EQ.1) THEN
            ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WPR
            PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WPR
          ELSEIF (ITYP.EQ.2) THEN
            ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WPR
            PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WPR
          ELSEIF (ITYP.EQ.3) THEN
            ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WPR
            PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WPR
          ENDIF
        ENDIF
C
90      CONTINUE
        IF (NADSI.GE.1) CALL UPSUSR (WPR,2)
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 4'
c slmod end
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
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN ? - 2',E0
c slmod end
      IF (ICOL.EQ.1) RETURN 3
C
C  ..........................................................................
C
C  NOW DEAL WITH REFLECTED AND/OR SPUTTERED PARTICLES.
C  ..........................................................................
C
      IF (ISRS(ISPZ,MSURF).NE.0..OR.ISRC(ISPZ,MSURF).NE.0.) THEN
C
C  FOLLOW BOTH: SPUTTERED AND REFLECTED PARTICLES
C  IN PRACTICE, FOLLOW ONLY ONE SECONDARY WITH WEIGHT ADJUSTMENT
C  AND EQUAL CHANCES FOR BOTH THE SPUTTERED AND REFLECTED PARTICLE
C  EXPECTED WEIGHT OF NEW TEST FLIGHT: WSPUT+WREFL
C
        WSPUT=WGHTSP+WGHTSS
        PSPUT=0.
        IF (WSPUT.GT.0.) PSPUT=0.5
      ELSE
C
C  DO NOT FOLLOW SPUTTERED PARTICLES
C
        WSPUT=0.
        PSPUT=0.
      ENDIF
C
C  REFLECTION FROM SURFACE
C
C  .......................................
C  .                                     .
C  .  MOLECULE REFLECTION MODEL 600--699 .
C  .......................................
C
C
600   CONTINUE
C
      IF (ITYP.EQ.2) THEN
C
C  TEST FOR REFLECTION
C
c slmod begin
        IF (ISRT(ISPZ,MSURF).EQ.0) THEN
          SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
          LGPART=.FALSE.
          RETURN
        ELSEIF (WEIGHT.GE.WMINS) THEN
c
c        IF (WEIGHT.GE.WMINS) THEN
c slmod end
C  WITH SUPPRESSION OF ABSORPTION
          WABS=WEIGHT*(1.D0-RECYCT(ISPZ,MSURF))
          IF (WABS.GT.0.D0) THEN
            SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WABS
            WEIGHT=WEIGHT-WABS
          ENDIF
          IF (WEIGHT.GT.EPS30) GOTO 610
          LGPART=.FALSE.
c slmod begin - debug - tr
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,*) '   ESCAPE: RETURN 5'
c slmod end
          RETURN
        ELSE
C  NO SUPPRESSION OF ABSORPTION
          ZVZ=RANF_EIRENE( )
          IF (ZVZ.LT.RECYCT(ISPZ,MSURF)) GOTO 610
C  ABSORB THIS PARTICLE
          SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
          LGPART=.FALSE.
c slmod begin - debug - tr
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,*) '   ESCAPE: RETURN 6'
c slmod end
          RETURN
        ENDIF
C
610     CONTINUE
C
C  NEW SPECIES: AGAIN MOLECULE
C
        IMOL=ISRT(ISPZ,MSURF)
        IF (IMOL.LT.1.OR.IMOL.GT.NMOLI) THEN
          FR2=RANF_EIRENE( )
          DO 621 I=1,NMOLI
            IMOL=I
            IF (FR2.LE.DMOL(IMOL)) GOTO 622
621       CONTINUE
          WRITE (6,*) 'SPECIES INDEX OUT OF RANGE IN ESCAPE '
          WRITE (6,*) 'IMOL, MSURF ',IMOL,MSURF
          CALL EXIT
622       CONTINUE
        ENDIF
C
        ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
C
c slmod begin - debug - tr


        E0TERM=EWALL(MSURF)

      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: A E0= ',E0,E0TERM


c        IF (0.667D0*E0.LT.-EWALL(MSURF)) THEN
c          E0TERM=EWALL(MSURF)
c        ELSE
c          E0TERM=0.95D0*E0
c        ENDIF

c        E0TERM=DMAX1(0.65D0*E0,-1.5D0*EWALL(MSURF))

c        WRITE(0,'(A,1P,3E10.2,0P)') 
c     .        'E?:',E0term,e0,-1.5D0*EWALL(MSURF)


c        IF (MSURF.LE.6.OR.MSURF.GT.NLIM) THEN
c          E0TERM=DMAX1(0.95D0*E0,-1.5D0*EWALL(MSURF))
c        ELSE
c          E0TERM=EWALL(MSURF)
c        ENDIF
        



c        E0TERM=E0

c
c        E0TERM=EWALL(MSURF)
c slmod end
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
          CALL VELOCS (TW,0.D0,0.D0,0.D0,0.D0,0.D0,RSQDVM(IMOL),
     .                 CVRSSM(IMOL),
     .                -CRTX,-CRTY,-CRTZ,
     .                 E0,VELX,VELY,VELZ,VEL)
          E0_MEAN=2.*TW
        ELSE
          WRITE (6,*) 'ERROR IN ESCAPE, EXIT CALLED '
          CALL EXIT
        ENDIF
        WREFL=WEIGHT
C
C  ...............................................
C  .                                             .
C  .  REFLECTION MODEL FOR ATOMS OR ATOMIC IONS  .
C  ...............................................
C
      ELSEIF (ITYP.EQ.1.OR.ITYP.EQ.3) THEN
C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: C E0= ',E0

     

        IF (MSURF.GE.3.AND.MSURF.LE.6) THEN
          HOLDE0=E0
        ENDIF
c slmod end
        CALL REFLC1 (WMINS,FMASS,FCHAR,NPRT(ISPZ),
     .               ISRF(ISPZ,MSURF),ISRT(ISPZ,MSURF))
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: D E0= ',E0
        IF (MSURF.GE.3.AND.MSURF.LE.6.AND.ityp.EQ.1) THEN 
          ndeltae0 = ndeltae0+1.0D0
          deltae0=deltae0+(e0/holde0)


c     .    WRITE(0,*) '   ESCAPE: D E0= ',ityp
        ENDIF
c slmod end
        ISPZ=ISPEZ(ITYP,IATM,IMOL,IION,IPLS)
        WREFL=WEIGHT
        IF (.NOT.LGPART) THEN
          WREFL=0.
        ENDIF
      ENDIF
C
C
C  DECIDE: CONTINUE WITH REFLECTED OR SPUTTERED PARTICLE
C
      ZEP1=RANF_EIRENE( )
      IF (ZEP1.LE.PSPUT) THEN
C  DECISION IS MADE: SPUTTERED PARTICLE WILL BE FOLLOWED NEXT
        IF (RANF_EIRENE( ).LE.WGHTSP/WSPUT) THEN
C  PHYSICALLY SPUTTERED PARTICLE
          ISPZ=ISSPTP
          ITYP=ISPEZI(ISPZ,0)
          IATM=ISPEZI(ISPZ,1)
          IMOL=ISPEZI(ISPZ,2)
          IION=ISPEZI(ISPZ,3)
          IPLS=ISPEZI(ISPZ,4)
          E0=ESPTP
          VEL=VSPTP
          VELX=VXSPTP
          VELY=VYSPTP
          VELZ=VZSPTP
C
          WEIGHT=WSPUT/PSPUT
          LGPART=NLSPUT
        ELSE
C  CHEMICALLY SPUTTERED PARTICLE
          ISPZ=ISSPTC
          ITYP=ISPEZI(ISPZ,0)
          IATM=ISPEZI(ISPZ,1)
          IMOL=ISPEZI(ISPZ,2)
          IION=ISPEZI(ISPZ,3)
          IPLS=ISPEZI(ISPZ,4)
          E0=ESPTC
          VEL=VSPTC
          VELX=VXSPTC
          VELY=VYSPTC
          VELZ=VZSPTC
C
          WEIGHT=WSPUT/PSPUT
          LGPART=NLSPUT
        ENDIF
      ELSE
C  DECISION IS MADE: REFLECTED PARTICLE WILL BE FOLLOWED NEXT
        WEIGHT=WREFL/(1.-PSPUT)
      ENDIF
C
C  UPDATE REFLECTED PARTICLE AND ENERGY FLUX
C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN ? - 3',E0
c slmod end
      IF (.NOT.LGPART) RETURN
C
C  ITOLD.NE.ITNEW=ITYP POSSIBLE
C
      IF (ITYP.EQ.1) THEN
        LOGATM(IATM,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WEIGHT
          ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WEIGHT
            ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 7'
c slmod end
          RETURN 1
        ELSEIF (ITOLD.EQ.2) THEN
          PRFMAT(IATM,MSURF)=PRFMAT(IATM,MSURF)+WEIGHT
          ERFMAT(IATM,MSURF)=ERFMAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            PRFMAT(IATM,MSURFG)=PRFMAT(IATM,MSURFG)+WEIGHT
            ERFMAT(IATM,MSURFG)=ERFMAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 8'
c slmod end
          RETURN 1
        ELSEIF (ITOLD.EQ.3) THEN
          PRFIAT(IATM,MSURF)=PRFIAT(IATM,MSURF)+WEIGHT
          ERFIAT(IATM,MSURF)=ERFIAT(IATM,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            PRFIAT(IATM,MSURFG)=PRFIAT(IATM,MSURFG)+WEIGHT
            ERFIAT(IATM,MSURFG)=ERFIAT(IATM,MSURFG)+E0*WEIGHT
          ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 9'
c slmod end
          RETURN
        ENDIF
      ELSEIF (ITYP.EQ.2) THEN
        LOGMOL(IMOL,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          PRFAML(IMOL,MSURF)=PRFAML(IMOL,MSURF)+WEIGHT
          ERFAML(IMOL,MSURF)=ERFAML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            PRFAML(IMOL,MSURFG)=PRFAML(IMOL,MSURFG)+WEIGHT
            ERFAML(IMOL,MSURFG)=ERFAML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 10'
c slmod end
          RETURN 1
        ELSEIF (ITOLD.EQ.2) THEN
          PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WEIGHT
          ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WEIGHT
            ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 11'
c slmod end
          RETURN 1
        ELSEIF (ITOLD.EQ.3) THEN
          PRFIML(IMOL,MSURF)=PRFIML(IMOL,MSURF)+WEIGHT
          ERFIML(IMOL,MSURF)=ERFIML(IMOL,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            PRFIML(IMOL,MSURFG)=PRFIML(IMOL,MSURFG)+WEIGHT
            ERFIML(IMOL,MSURFG)=ERFIML(IMOL,MSURFG)+E0*WEIGHT
          ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 12'
c slmod end
          RETURN
        ENDIF
      ELSEIF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          PRFAIO(IION,MSURF)=PRFAIO(IION,MSURF)+WEIGHT
          ERFAIO(IION,MSURF)=ERFAIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            PRFAIO(IION,MSURFG)=PRFAIO(IION,MSURFG)+WEIGHT
            ERFAIO(IION,MSURFG)=ERFAIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 13'
c slmod end
          RETURN
        ELSEIF (ITOLD.EQ.2) THEN
          PRFMIO(IION,MSURF)=PRFMIO(IION,MSURF)+WEIGHT
          ERFMIO(IION,MSURF)=ERFMIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            PRFMIO(IION,MSURFG)=PRFMIO(IION,MSURFG)+WEIGHT
            ERFMIO(IION,MSURFG)=ERFMIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 14'
c slmod end
          RETURN
        ELSEIF (ITOLD.EQ.3) THEN
          PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WEIGHT
          ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WEIGHT
          IF (MSURFG.GT.0) THEN
            PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WEIGHT
            ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WEIGHT
          ENDIF
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 15'
c slmod end
          RETURN 1
        ENDIF
      ENDIF
      IF (NADSI.GE.1) CALL UPSUSR (WPR,2)
C
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,*) '   ESCAPE: RETURN 16'
c slmod end
      RETURN
      END
C
C
      SUBROUTINE REFLEC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  REFLECT ESCAPING ATOMS OR IONS
C  INPUT:
C       ILREF = 1  DATABASE REFLECTION MODEL, W.ECKSTEIN, D.B.HEIFETZ,
C                  IPP 9/59 (1986)
C       ILREF = 2  MODIFIED BEHRISCH MATRIX, R. BEHRISCH, ERICE SUMMER
C                  SCHOOL 1976
C       ILREF = 3  USER SUPPLIED REFLECTION MODEL, CALL: RF1USR
C
C       ITYP  = 1  INCIDENT ATOM
C       ITYP  = 3  INCIDENT TEST ION
C       ITYP  = 4  INCIDENT BULK ION
C  OUTPUT:
C     LGPART= TRUE AND:
C       ITYP = 1  ATOM IATM IS RETURNED TO CALLING PROGRAM
C       ITYP = 2  MOLECULE IMOL IS RETURNED TO CALLING PROGRAM
C       ITYP = 3  TEST ION  IION IS RETURNED TO CALLING PROGRAM
C     LGPART= FALSE  NO PARTICLE IS RETURNED (ABSORBTION)
C       ITYP = 0
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CZT1'
      INCLUDE 'CLOGAU'
      INCLUDE 'COMUSR'
      INCLUDE 'CCONA'
      INCLUDE 'CRAND'
      INCLUDE 'CTRCEI'
      INCLUDE 'CREF'
      INCLUDE 'CLGIN'
      INCLUDE 'CADGEO'
      INCLUDE 'CESTIM'
      INCLUDE 'CSPEI'
      DIMENSION HFTR0(NHD1,NHD2,NHD6),
     .          HFTR1(NHD1,NHD2,NHD3,NHD6),
     .          HFTR2(NHD1,NHD2,NHD3,NHD4,NHD6),
     .          HFTR3(NHD1,NHD2,NHD3,NHD4,NHD5,NHD6),
     .          HFTR3F(NHD5)
      EQUIVALENCE (RWK(NID2+1),HFTR0(1,1,1)),
     .            (RWK(NID2+1+NH0),HFTR1(1,1,1,1)),
     .            (RWK(NID2+1+NH0+NH1),HFTR2(1,1,1,1,1)),
     .            (RWK(NID2+1+NH0+NH1+NH2),HFTR3(1,1,1,1,1,1))
C---------------------------------------------------------------------
C
      DIMENSION
     .  ZRANGE(0:12),ZDE(12),ZDEL(12),ZENGY(0:12),ZR(0:12),ZIDE(12,12),
     .  ZIDED(12,12),XSP(12),YSP(12),ASP(12),BSP(12),CSP(12),DSP(12),
     .  E0AV(0:12),QUOTR(0:11),QUOTE(0:11),ERDC(NHD6)
      LOGICAL NLDATA,NLBEHR
      SAVE
C  SIZE OF "BEHRISCH TABLES"
      DATA IDIM/12/
C  ENERGY , ABSZISSA FOR REFLECTION PROBABILITY, H INCIDENT ON FE
      DATA
     .   ZENGY/0.,4.64,10.0,21.5,46.4,100.0,215.4,464.1,1000.0,2154.3,
     .         4641.3,10000.0,21543.0/
C  REFLECTION PROBABILITY RPROB(ENERGY)= ZR(ZENGY)
      DATA
     .   ZR/1.,0.9,0.8,0.7,0.62,0.543,0.46,0.37,0.29,0.21,0.14,
     .      0.095,0.04/
C  ENERGY RANGE FOR ENERGY DISTRIBUTION, LAST CELL IS: ZENGY
C  I.E. ABSZISSA FOR ENERGY-DISTRIBUTION-FUNCTIONS, H INCIDENT ON FE
      DATA
     .   ZRANGE/0.,6.81,14.7,31.63,68.1,146.8,316.3,681.9,1468.0,
     .          3162.0,6813.0,14678.0,31630./
C  DISTRIBUTION FUNCTIONS ZIDE(ZRANGE) , ONE FOR EACH ZENGY
      DATA ZIDE
     . /12*1.,
     .  0.2,11*1.,
     .  0.1,0.2,10*1.,
     .  0.025,0.05,0.35,9*1.,
     .  0.012,0.025,0.1,0.45,8*1.,
     .  0.01,0.02,0.05,0.15,0.55,7*1.,
     .  0.002,0.005,0.02,0.055,0.175,0.625,6*1.,
     .  0.001,0.003,0.011,0.029,0.079,0.224,0.699,5*1.,
     .  0.001,0.003,0.007,0.016,0.04,0.105,0.301,0.771,4*1.,
     .  0.,0.001,0.003,0.007,0.018,0.050,0.14,0.40,0.83,3*1.,
     .  0.,0.001,0.003,0.006,0.011,0.025,0.073,0.215,0.505,0.865,2*1.,
     .  0.,0.001,0.003,0.005,0.009,0.015,0.035,0.105,0.305,0.6,0.9,1./
C---------------------------------------------------------------------
      DATA CON/0.4685/,EOQ/14.39/,ZWDR/0.666667/,IFIRST/0/,ICOUNT/0/
      EREDC(XMTT,XCTT,XMPP,XCPP)=CON/EOQ*XMTT/((XMPP+XMTT)*XCPP*XCTT*
     .                       SQRT(XCPP**ZWDR+XCTT**ZWDR))
C
C  INITIALIZE SURFACE REFLECTION MODELS
C
      ENTRY REFLC0
C
      IF (IFIRST.EQ.1) RETURN
      IFIRST=1
C
      NLDATA=.FALSE.
      NLBEHR=.FALSE.
      DO 1 J=1,NLIMPS
        NLDATA=NLDATA.OR.(ILREF(J).EQ.1)
        NLBEHR=NLBEHR.OR.(ILREF(J).EQ.2)
1     CONTINUE
C
C
C  SET ADDITIONAL DATA FOR "DATABASE REFLECTION MODEL"
C
      IF (NLDATA) THEN
C
        IF (NLTRIM) THEN
c          IF (output)
c     .    WRITE(6,*) 'MARK: REFLC0: LTRMOL= ',LTRMOL
          IF (LTRMOL) THEN
            CALL REFDAT(TM,TC,WM,WC)
          ELSE
            CALL RDTRIM
          ENDIF
        ELSE
          WRITE (6,*) 'INPUT ERROR FOR LOCAL REFLECTION MODEL'
          WRITE (6,*) 'DATABASE REFLECTION MODEL REQUIRED BUT '
          WRITE (6,*) 'NLTRIM IS NOT SET TRUE'
          CALL EXIT
        ENDIF
C
C  SET FACTORS FOR REDUCED ENERGY SCALING FOR ALL TARGET/PROJECTILE
C  COMBINATIONS AVAILABLE IN DATABASE MODEL
        DO 3 J=1,NFLR
          ERDC(J)=EREDC(WM(J),WC(J),TM(J),TC(J))
3       CONTINUE
C  SET UNIFORM DISTRIBUTION OF AZIMUTAL ANGLE FOR DATABASE MODEL
C  FOR PERPENDICULAR INCIDENCE (INDW=1)
        DO 4 INDR3=1,INR
          HFTR3F(INDR3)=COS(PIA*(1.-RAAR(INDR3)))
4       CONTINUE
C
        IF (TRCREF) THEN
          DO 5 J=1,NFLR
            CALL LEER(1)
            WRITE (6,*) 'DATABASE REFLECTION MODEL DEFINED FOR:'
            WRITE (6,*) 'IFILE =                      ',J
            WRITE (6,*) 'TARGET MASS NUMBER =         ',WM(J)
            WRITE (6,*) 'TARGET NUCL. CHARGE NUMBER = ',WC(J)
            WRITE (6,*) 'PROJECTIL MASS NUMBER =      ',TM(J)
            WRITE (6,*) 'PRJTL. NUCL. CHARGE NUMBER = ',TC(J)
            WRITE (6,*) 'REDUCED ENERGY FACTOR ERDC = ',ERDC(J)
5         CONTINUE
          CALL LEER(2)
        ENDIF
      ELSE
        INE=1
        INW=1
        INR=1
        NFLR=1
      ENDIF
C
C  SET ADDITIONAL DATA FOR "BEHRISCH MATRIX REFLECTION MODEL"
C
      IF (NLBEHR) THEN
C
C  "BEHRISCH MODEL" REFLECTION PROBABILITY ZR(ZENGY)
C        FOR ENERGIES BELOW ERCUT CAN BE MODIFIED;
C        RPROB0 IS THE NEW (HYPOTHETICAL) REFL. PROB AT ZENGY(0)=0. (EV)
C        (THERE IS NO REFLECTION MODEL BUT ONLY A THERMAL PARTICLE
C        MODEL CALLED FOR E0 BELOW ERMIN)
C
C  NO MODIFICATION FOR ERCUT.LT.0.
C  ZR(0) -- ZR(NRE) ARE MODIFIED FOR ERCUT.GE.ZENGY(0)
        IF (ERCUT.GT.ZENGY(0).AND.ERCUT.LT.ZENGY(1)) THEN
          ZR(0)=RPROB0
        ELSEIF (ERCUT.GE.ZENGY(1)) THEN
          NRI=2
          NRE=LEARCA(ERCUT,ZENGY,1,13,1,'REFLEC (1)  ')-1
          NREP=NRE+1
          DX=ERCUT-ZENGY(NRE)
          XSP(1)=ZENGY(0)
          YSP(1)=RPROB0
          XSP(2)=ERCUT
          YSP(2)=(ZR(NREP)-ZR(NRE))/(ZENGY(NREP)-ZENGY(NRE))*DX+ZR(NRE)
C  SET SPLINE DATA FOR NEW REFLECTION PROBABLITY ZR(ZENGY)
          DO 6 J=NREP,12
            NRI=NRI+1
            XSP(NRI)=ZENGY(J)
            YSP(NRI)=ZR(J)
6         CONTINUE
          CALL SPLINE(XSP,YSP,NRI,ASP,BSP,CSP,DSP)
C   SET THE NEW ZR ON GRID ZENGY(J) FROM ZENGY(0) TO ZENGY(NRE)
          ZR(0)=RPROB0
          DO 7 J=0,NRE
            DX=ZENGY(J)-XSP(1)
            ZR(J)=((DSP(1)*DX+CSP(1))*DX+BSP(1))*DX+ASP(1)
7         CONTINUE
        ENDIF
C
C
C  THE BEHRISCH-REFLECTION DATA ARE GIVEN FOR H INCIDENT ON FE
C  CONVERT TO REDUCED ENERGY
C
C  CHARGE NUMBERS: STAINLESS STEEL
        XCFE=26.
C  MASS NUMBERS  : STAINLESS STEEL
        XMFE=56.
C  CHARGE NUMBER: HYDROGEN
        XCH=1.
C  MASS NUMBER  : HYDROGEN
        XMH=1.
C
        EPSHFE=EREDC(XMFE,XCFE,XMH,XCH)
C
        DO 12 J=1,12
          ZRANGE(J)=ZRANGE(J)*EPSHFE
          ZDE(J)=ZRANGE(J)-ZRANGE(J-1)
          ZENGY(J)=ZENGY(J)*EPSHFE
12      CONTINUE
        DO 13 I=1,12
          ZIDED(1,I)=ZIDE(1,I)
          ZDEL(I)=ZENGY(I)-ZRANGE(I-1)
          DO 14 J=2,12
            ZIDED(J,I)=ZIDE(J,I)-ZIDE(J-1,I)
14        CONTINUE
13      CONTINUE
        ZDEL(1)=0.
C
C  SET MEAN ENERGY OF REFLECTED PARTICLES FROM STOCHASTIC MATRIX
C
        E0AV(0)=0.
        DO 15 J=1,12
C  LAST BOX
          E0AV(J)=ZIDED(J,J)*(ZENGY(J)-ZDEL(J)/2.)
C   OTHER BOXES
          DO 16 I=1,J-1
16          E0AV(J)=E0AV(J)+ZIDED(I,J)*(ZRANGE(I)-ZDE(I)/2.)
15      CONTINUE
C
C  SET SOME CONSTANTS TO SPEED UP LINEAR INTERPOLATION IN
C  BEHRISCH REFLECTION DATA
C
        DO 17 J=0,11
          JP=J+1
          QUOTR(J)=(ZR(JP)-ZR(J))/(ZENGY(JP)-ZENGY(J))
          QUOTE(J)=(E0AV(JP)-E0AV(J))/(ZENGY(JP)-ZENGY(J))
17      CONTINUE
C
        IF (TRCREF) THEN
          CALL LEER(1)
          WRITE (6,*) 'REFLECTION DATA FROM BEHRISCH-MATRIX'
          WRITE (6,*) 'IMP. ENERGY (RED), REF. PROB, MEAN REFL. ENERGY'
          DO 19 J=0,12
            CALL MASR3('                        ',ZRANGE(J),ZR(J),
     .                                            E0AV(J))
19        CONTINUE
          CALL LEER(2)
        ENDIF
C
      ENDIF
C
C  PRINTOUT REFLECTION PROPERTIES OF SURFACES
C
      IF (TRCREF) THEN
        WRITE (6,*) 'ADDITIONAL SURFACES, THAT ARE NOT 100% RECYCLING'
        WRITE (6,*) 'FOR ALL SPECIES'
        ICOUNT=0
        DO 20 ILIM=1,NLIMI
          DO 21 ISP=1,NSPZ
            IF (RECYCT(ISP,ILIM).LT.1.D0) THEN
              ICOUNT=ICOUNT+1
              IF (ICOUNT.EQ.1) WRITE (6,*) 'ISPZ,ILIM,RECYCT'
              WRITE (6,*) ISP,ILIM,RECYCT(ISP,ILIM)
            ENDIF
21        CONTINUE
20      CONTINUE
        IF (ICOUNT.EQ.0) WRITE (6,*) 'NONE '
        WRITE (6,*) 'STANDARD SURFACES, THAT ARE NOT 100% REFLECTING'
        WRITE (6,*) 'FOR ALL SPECIES'
        ICOUNT=0
        DO 22 ISTS=1,NSTSI
          DO 23 ISP=1,NSPZ
            IF (RECYCT(ISP,NLIM+ISTS).LT.1.D0) THEN
              ICOUNT=ICOUNT+1
              IF (ICOUNT.EQ.1) WRITE (6,*) 'ISPZ,ISTS,RECYCT'
              WRITE (6,*) ISP,ISTS,RECYCT(ISP,NLIM+ISTS),NLIM,
     .                    ILIIN(NLIM+ISTS)
            ENDIF
23        CONTINUE
22      CONTINUE
        IF (ICOUNT.EQ.0) WRITE (6,*) 'NONE '
        ICOUNT=0
        CALL LEER(2)
      ENDIF
C
      CALL RF0USR
C
      RETURN
C
      ENTRY REFLC1 (WMIN,XMP,XCP,NPRIN,IGASF,IGAST)
C
C  SURFACE NUMBER  : MSURF (MSURF=0: DEFAULT MODEL)
C  SPECIES INDEX   : ISPZ
C
      MODREF=ILREF(MSURF)
      XMW=ZNML(MSURF)
      XCW=ZNCL(MSURF)
      E0TERM=EWALL(MSURF)
      EBIND=EWBIN(MSURF)
      PRFCF=RECYCF(ISPZ,MSURF)
      PRFCT=RECYCT(ISPZ,MSURF)
      EXPP=EXPPL(ISPZ,MSURF)
      EXPE=EXPEL(ISPZ,MSURF)
      EXPI=EXPIL(ISPZ,MSURF)
      RINTG=RINTEG(MSURF)
      EINTG=EINTEG(MSURF)
      AINTG=AINTEG(MSURF)
      ISPZO=ISPZ
c      IF (output)
c     .WRITE(6,*) 'MARK: REFLC1: MSURF,MODREF= ',MSURF,MODREF
C
C  SET EMIN AND EMAX FOR ENERGY SAMPLING FROM DATABASE (MODREF=1)
C
      EMINR=E0TERM
      IF (E0TERM.LT.0.D0) EMINR=-2.*E0TERM
      EMAXR=E0
C
C   TENTATIVELY ASSUME  REFLECTION
      LGPART=.TRUE.
C   COSINE OF ANGLE OF INCIDENCE
      COSIN=VELX*CRTX+VELY*CRTY+VELZ*CRTZ
      IF (COSIN.LT.0.D0) GOTO 990
C
C   NO REFLECTION OF FAST ATOMS FOR INCIDENT ENERGY BELOW ERMIN
C                               OR IF IGASF=0
C
c      IF (ityp.EQ.4) WRITE(0,*) 'debug e0:',e0,ermin,igasf
      IF (E0.LE.ERMIN.OR.IGASF.EQ.0) THEN
C
C   THERMAL PARTICLE MODEL IS CALLED
C
        RPROB=0.
        WFAC=0.
        COSIN=1.
C  RELATIVE FRACTION OF COSINE VS. SPECULAR REFLECTION
        F1=1.
        F2=0.
        FR1=RANF_EIRENE( )
C       IF (FR1.GE.RPROB) THEN
          IF (IGAST) 500,700,600
C       ENDIF
      ENDIF
C
C   FACTOR FOR CONVERSION TO REDUCED ENERGY
      EREDUC=EREDC(XMW,XCW,XMP,XCP)
C
C
C   MODREF=1: "DATABASE REFLECTION MODEL" (TRIM)
C   MODREF=2: "BEHRISCH-MATRIX"
C   MODREF=3: "USER SUPPLIED REFLECTION MODEL"
C
      IF (MODREF.EQ.1) THEN
        GOTO 100
      ELSEIF (MODREF.EQ.2) THEN
        GOTO 200
      ELSEIF (MODREF.GE.3) THEN
        CALL RF1USR (XMW,XCW,XMP,XCP,IGASF,IGAST,F1,F2,EXPI,
     .               RPROB,E0TERM,*400,*500,*600,*700)
        RETURN
      ENDIF
C
C  DATABASE REFLECTION MODEL STARTS HERE
C
100   CONTINUE
C
C   CHECK IF WALL REFLECTION DATA FOR IATM/IION INCIDENT ON
C   XWALL/ZWALL ARE AVAILABLE
C
      EQTO=1.D40
      EFCT=1.
      DO 120 IFILE=1,NFLR
        IF (output.AND.trcref)
     .  WRITE(6,'(A,2I4,1P,2E15.7,0P)') ' MARK: REFLC1: IFILE ',
     .    IFILE,NFLR,ERDC(IFILE),EREDUC
        IF (ABS(ERDC(IFILE)-EREDUC).LE.EPS12) GOTO 130
        EQT=EREDUC/ERDC(IFILE)
        ETEST=ABS(EQT-1.)
        IF (ETEST.LT.EQTO) THEN
          ISAVE=IFILE
          EQTO=ETEST
          EQSAVE=EQT
          IF (output.AND.trcref)
     .    WRITE(6,'(A,1P,2E12.4,0P,I4)') ' MARK: REFLC1: EQTO ',
     .      ETEST,EQTO,ISAVE
        ENDIF
120   CONTINUE
c      IF (output)
c     .WRITE(6,*) 'MARK: REFLC1: IFILE= ',IFILE
      IF (ICOUNT.LT.5.AND.TRCREF) THEN
        WRITE (6,*) 'TRIM-REFLECTION DATA REQUIRED BUT NOT'
        WRITE (6,*) 'AVAILABLE FOR THE TARGET-PROJECTIL SYSTEM:'
        WRITE (6,*) 'XMWALL,XCWALL,XMPART,XCPART '
        WRITE (6,*)  XMW,XCW,XMP,XCP
        WRITE (6,*) 'EREDUC = ',EREDUC
        WRITE (6,*) 'THE REDUCED ENERGY FORMULAS ARE APPLIED WITH'
        WRITE (6,*) 'THE DATA FOR THE TARGET-PROJECTIL SYSTEM:'
        WRITE (6,*) 'J,WM(J),WC(J),TM(J),TC(J) '
        WRITE (6,*)  ISAVE,WM(ISAVE),WC(ISAVE),TM(ISAVE),TC(ISAVE)
        WRITE (6,*) 'ERDC(J),J=1,NFLR = ',(ERDC(J),J=1,NFLR)
        CALL LEER(1)
        ICOUNT=ICOUNT+1
      ENDIF
      IFILE=ISAVE
      EFCT=EQSAVE
C
      E0=E0*EFCT
C
130   CONTINUE
C
C  FIND INDICES FOR INCIDENT ENERGY AND ANGLE: INDE, INDW, R01, R02
C
      DO 102 I=2,INEM
        INDEP=I
        IF (E0.LE.ENAR(I)) GOTO 101
102   CONTINUE
      INDEP=INE
101   INDE=INDEP-1
C
      DO 103 I=2,INWM
         INDWP=I
         IF (COSIN.GE.WIAR(I)) GOTO 104
103      CONTINUE
      INDWP=INW
104   INDW=INDWP-1
C
c      IF (output) THEN
c        WRITE(6,*) 'MARK: REFCT1: RO1: ',INDE,E0,ENAR(INDE),DENAR(INDE)
c        WRITE(6,*) 'MARK: REFCT1: RO2: ',INDW,COSIN,WIAR(INDW),
c     .                                   DWIAR(INDW)
c        WRITE(6,*) 'MARK: REFCT1: HFTR0: ',
c     .              INDE,INDW,IFILE,HFTR0(INDE,INDW,IFILE)
c        WRITE(6,*) 'MARK: REFCT1: HFTR0: ',
c     .              INDE,INDWP,IFILE,HFTR0(INDE,INDWP,IFILE)
c      ENDIF
      RO1=(E0-ENAR(INDE))*DENAR(INDE)
      RO2=(COSIN-WIAR(INDW))*DWIAR(INDW)
C
C  REFLECTION PROBALITY: RPROB
C
      IF (RINTG.GT.0.D0) THEN
        RPROB=MIN(PRFCT,RINTG)
      ELSE
        RF1=HFTR0(INDE,INDW,IFILE)
        RF1=RF1+RO1*(HFTR0(INDEP,INDW,IFILE)-RF1)
        RF2=HFTR0(INDE,INDWP,IFILE)
        RF2=RF2+RO1*(HFTR0(INDEP,INDWP,IFILE)-RF2)
C
        RPROB=RF1+RO2*(RF2-RF1)
c slmod begin - tr
c        IF (PRFCF.LT.0.0) THEN
c          RPROB=MIN(-PRFCF,PRFCT)
c        ELSE
c          RPROB=MIN(RPROB*PRFCF,PRFCT)
c        ENDIF
c
        RPROB=MIN(RPROB*PRFCF,PRFCT)
c slmod end
        IF (output) THEN
          WRITE(6,'(A,2I4,F10.4,I4)')
     .      ' MARK: REFCT1: RF1: ',INDE,INDW,RF1,IFILE
          WRITE(6,'(A,I4,F10.4,I4)')
     .      ' MARK: REFCT1: RF2: ',INDEP,RF2,IFILE
        ENDIF
      ENDIF
C
C   DECIDE IF PARTICLE IS TO BE REFLECTED OR IF THE "THERMAL
C   PARTICLE-MODEL" IS CALLED
C
      WFAC=1.
      FR1=RANF_EIRENE( )
      IATM=IGASF
C  THERMAL PARTICLE MODEL
      IF (output)
     .  WRITE(6,'(A,I6,2F14.6)') 
     .   ' MARK: REFLC1: IGAST,FR1,RPROB= ',IGAST,FR1,RPROB
     
c...TEMP
      IF (opttest.EQ.1) fr1 = -1.0D0

c      IF (FR1.GE.RPROB) 

c        WRITE(0,'(A,2F12.4,2I6,L3,I6,F12.4)') 
c     .   'MARK:RPROB=',rprob,fr1,msurf,
c     .   ityp,FR1.GE.RPROB,IGAST,e0

c      IF (ityp.EQ.4) THEN
c        icnt = icnt + 1
c        if (fr1.LT.rprob) ifast = ifast + 1 
c        if (fr1.GT.rprob) islow = islow + 1 
c        IF (mod(icnt,100).EQ.0) THEN
c          WRITE(0,*) 'debug rprob:',e0
c          WRITE(0,*) 'debug rprob:',fr1,rprob
c          WRITE(0,*) 'debug rprob:',ifast,islow,REAL(ifast)/REAL(islow)
c        endif
c      ENDIF
      IF (FR1.GE.RPROB) THEN
        IF (IGAST) 500,700,600
      ENDIF
C
C  ENERGY OF REFLECTED PARTICLE
C
      ZEP1=RANF_EIRENE( )
      DO 105 I=2,INRM
        INDR1P=I
        IF (ZEP1.LE.RAAR(I)) GOTO 106
105   CONTINUE
      INDR1P=INR
106   INDR1=INDR1P-1
C
      RO3=(ZEP1-RAAR(INDR1))*DRAAR(INDR1)
C
      IF (EINTG.GT.0.D0) THEN
C  CONSTANT ENERGY REFLECTION COEFFICIENT
        E0=E0*EINTG
C     ELSEIF (EINTG.LT.0.D0) THEN
C  E0 FROM MEAN ENERGY MODEL
      ELSE
C  E0 FROM STOCHASTIC MATRIX
        RF1=HFTR1(INDE,INDW,INDR1,IFILE)
        RF1=RF1+RO1*(HFTR1(INDEP,INDW,INDR1,IFILE)-RF1)
        RF2=HFTR1(INDE,INDWP,INDR1,IFILE)
        RF2=RF2+RO1*(HFTR1(INDEP,INDWP,INDR1,IFILE)-RF2)
        RF3=HFTR1(INDE,INDW,INDR1P,IFILE)
        RF3=RF3+RO1*(HFTR1(INDEP,INDW,INDR1P,IFILE)-RF3)
        RF4=HFTR1(INDE,INDWP,INDR1P,IFILE)
        RF4=RF4+RO1*(HFTR1(INDEP,INDWP,INDR1P,IFILE)-RF4)
C
        RFF1=RF1+RO2*(RF2-RF1)
        RFF2=RF3+RO2*(RF4-RF3)
C
        E0=RFF1+RO3*(RFF2-RFF1)
        E0=MAX(E0,EMINR)
        E0=MIN(E0,EMAXR)
      ENDIF
C
      E0=E0/EFCT
      VEL=RSQDVA(IATM)*SQRT(E0)
      ITYP=1
C
C  POLAR ANGLE OF REFLECTION
C
      IF (EXPI.EQ.0..OR.EXPI.GE.100.D0) THEN
        F1=1.
        F2=0.
        GOTO 400
      ENDIF
C
      ZEP1=RANF_EIRENE( )
      DO 107 I=2,INRM
        INDR2P=I
        IF (ZEP1.LE.RAAR(I)) GOTO 108
107   CONTINUE
      INDR2P=INR
108   INDR2=INDR2P-1
C
      RO4=(ZEP1-RAAR(INDR2))*DRAAR(INDR2)
C
      RF1=HFTR2(INDE,INDW,INDR1,INDR2,IFILE)
      RF1=RF1+RO1*(HFTR2(INDEP,INDW,INDR1,INDR2,IFILE)-RF1)
      RF2=HFTR2(INDE,INDWP,INDR1,INDR2,IFILE)
      RF2=RF2+RO1*(HFTR2(INDEP,INDWP,INDR1,INDR2,IFILE)-RF2)
      RF3=HFTR2(INDE,INDW,INDR1P,INDR2,IFILE)
      RF3=RF3+RO1*(HFTR2(INDEP,INDW,INDR1P,INDR2,IFILE)-RF3)
      RF4=HFTR2(INDE,INDWP,INDR1P,INDR2,IFILE)
      RF4=RF4+RO1*(HFTR2(INDEP,INDWP,INDR1P,INDR2,IFILE)-RF4)
      RF5=HFTR2(INDE,INDW,INDR1,INDR2P,IFILE)
      RF5=RF5+RO1*(HFTR2(INDEP,INDW,INDR1,INDR2P,IFILE)-RF5)
      RF6=HFTR2(INDE,INDWP,INDR1,INDR2P,IFILE)
      RF6=RF6+RO1*(HFTR2(INDEP,INDWP,INDR1,INDR2P,IFILE)-RF6)
      RF7=HFTR2(INDE,INDW,INDR1P,INDR2P,IFILE)
      RF7=RF7+RO1*(HFTR2(INDEP,INDW,INDR1P,INDR2P,IFILE)-RF7)
      RF8=HFTR2(INDE,INDWP,INDR1P,INDR2P,IFILE)
      RF8=RF8+RO1*(HFTR2(INDEP,INDWP,INDR1P,INDR2P,IFILE)-RF8)
C
      RFF1=RF1+RO2*(RF2-RF1)
      RFF2=RF3+RO2*(RF4-RF3)
      RFF3=RF5+RO2*(RF6-RF5)
      RFF4=RF7+RO2*(RF8-RF7)
C
      RFFF1=RFF1+RO3*(RFF2-RFF1)
      RFFF2=RFF3+RO3*(RFF4-RFF3)
C
      ZCPHI=RFFF1+RO4*(RFFF2-RFFF1)
C  LIMIT COSINE OF POLAR ANGLE TO 85. DEGREES
C  (I.E., 5 DEGREES AGAINST SURFACE TANGENTIAL PLANE)
      ZCPHI=MIN(0.999999D0,MAX(0.08716D0,ZCPHI))
      ZSPHI=SQRT(1.-ZCPHI*ZCPHI)
C
C  AZIMUTAL ANGLE OF REFLECTION
C
      ZEP1=RANF_EIRENE( )
      DO 109 I=2,INRM
         INDR3P=I
         IF (ZEP1.LE.RAAR(I)) GOTO 110
109      CONTINUE
      INDR3P=INR
110   INDR3=INDR3P-1
C
      RO5=(ZEP1-RAAR(INDR3))*DRAAR(INDR3)
C
      IF (INDW.EQ.1) THEN
        RF1=HFTR3F(INDR3)
        RF3=RF1
        RF5=RF1
        RF7=RF1
        RF9=HFTR3F(INDR3P)
        RF11=RF9
        RF13=RF9
        RF15=RF9
      ELSE
        RF1=HFTR3(INDE,INDW,INDR1,INDR2,INDR3,IFILE)
        RF1=RF1+RO1*(HFTR3(INDEP,INDW,INDR1,INDR2,INDR3,IFILE)-RF1)
        RF3=HFTR3(INDE,INDW,INDR1P,INDR2,INDR3,IFILE)
        RF3=RF3+RO1*(HFTR3(INDEP,INDW,INDR1P,INDR2,INDR3,IFILE)-RF3)
        RF5=HFTR3(INDE,INDW,INDR1,INDR2P,INDR3,IFILE)
        RF5=RF5+RO1*(HFTR3(INDEP,INDW,INDR1,INDR2P,INDR3,IFILE)-RF5)
        RF7=HFTR3(INDE,INDW,INDR1P,INDR2P,INDR3,IFILE)
        RF7=RF7+RO1*(HFTR3(INDEP,INDW,INDR1P,INDR2P,INDR3,IFILE)-RF7)
        RF9=HFTR3(INDE,INDW,INDR1,INDR2,INDR3P,IFILE)
        RF9=RF9+RO1*(HFTR3(INDEP,INDW,INDR1,INDR2,INDR3P,IFILE)-RF9)
        RF=HFTR3(INDE,INDW,INDR1P,INDR2,INDR3P,IFILE)
        RF11=RF+RO1*(HFTR3(INDEP,INDW,INDR1P,INDR2,INDR3P,IFILE)-RF)
        RF=HFTR3(INDE,INDW,INDR1,INDR2P,INDR3P,IFILE)
        RF13=RF+RO1*(HFTR3(INDEP,INDW,INDR1,INDR2P,INDR3P,IFILE)-RF)
        RF=HFTR3(INDE,INDW,INDR1P,INDR2P,INDR3P,IFILE)
        RF15=RF+RO1*(HFTR3(INDEP,INDW,INDR1P,INDR2P,INDR3P,IFILE)-RF)
      ENDIF
C
      RF2=HFTR3(INDE,INDWP,INDR1,INDR2,INDR3,IFILE)
      RF2=RF2+RO1*(HFTR3(INDEP,INDWP,INDR1,INDR2,INDR3,IFILE)-RF2)
      RF4=HFTR3(INDE,INDWP,INDR1P,INDR2,INDR3,IFILE)
      RF4=RF4+RO1*(HFTR3(INDEP,INDWP,INDR1P,INDR2,INDR3,IFILE)-RF4)
      RF6=HFTR3(INDE,INDWP,INDR1,INDR2P,INDR3,IFILE)
      RF6=RF6+RO1*(HFTR3(INDEP,INDWP,INDR1,INDR2P,INDR3,IFILE)-RF6)
      RF8=HFTR3(INDE,INDWP,INDR1P,INDR2P,INDR3,IFILE)
      RF8=RF8+RO1*(HFTR3(INDEP,INDWP,INDR1P,INDR2P,INDR3,IFILE)-RF8)
      RF10=HFTR3(INDE,INDWP,INDR1,INDR2,INDR3P,IFILE)
      RF10=RF10+RO1*(HFTR3(INDEP,INDWP,INDR1,INDR2,INDR3P,IFILE)-RF10)
      RF12=HFTR3(INDE,INDWP,INDR1P,INDR2,INDR3P,IFILE)
      RF12=RF12+RO1*(HFTR3(INDEP,INDWP,INDR1P,INDR2,INDR3P,IFILE)-RF12)
      RF14=HFTR3(INDE,INDWP,INDR1,INDR2P,INDR3P,IFILE)
      RF14=RF14+RO1*(HFTR3(INDEP,INDWP,INDR1,INDR2P,INDR3P,IFILE)-RF14)
      RF16=HFTR3(INDE,INDWP,INDR1P,INDR2P,INDR3P,IFILE)
      RF16=RF16+RO1*(HFTR3(INDEP,INDWP,INDR1P,INDR2P,INDR3P,IFILE)-RF16)
C
      RFF1=RF1+RO2*(RF2-RF1)
      RFF2=RF3+RO2*(RF4-RF3)
      RFF3=RF5+RO2*(RF6-RF5)
      RFF4=RF7+RO2*(RF8-RF7)
      RFF5=RF9+RO2*(RF10-RF9)
      RFF6=RF11+RO2*(RF12-RF11)
      RFF7=RF13+RO2*(RF14-RF13)
      RFF8=RF15+RO2*(RF16-RF15)
C
      RFFF1=RFF1+RO3*(RFF2-RFF1)
      RFFF2=RFF3+RO3*(RFF4-RFF3)
      RFFF3=RFF5+RO3*(RFF6-RFF5)
      RFFF4=RFF7+RO3*(RFF8-RFF7)
C
      RFFFF1=RFFF1+RO4*(RFFF2-RFFF1)
      RFFFF2=RFFF3+RO4*(RFFF4-RFFF3)
C
      ZCTHET=RFFFF1+RO5*(RFFFF2-RFFFF1)
      ZCTHET=MAX(-.999999D0,MIN(0.999999D0,ZCTHET))
      ZSTHET=SQRT(1.-ZCTHET*ZCTHET)
      ZSTHET=ZSTHET*SIGN(1.D0,(RANF_EIRENE( )-0.5))
C
      VX=-ZCPHI
      VY=ZSPHI*ZSTHET
      VZ=ZSPHI*ZCTHET
      IF (COSIN.GT.0.999999) THEN
        CALL ROTATF (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ)
      ELSE
        CALL ROTATE (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ,COSIN)
      ENDIF
      ITYP=1
      RETURN
C
C  MODIFIED BEHRISCH MATRIX MODEL STARTS HERE
C
200   CONTINUE
C
      E0=E0*EREDUC
C
C  DETERMINE INTERVAL FOR INCIDENT ENERGY: IRM, ED
C
      DO 201 J=1,IDIM
        IRANGE=J
        IF (ZENGY(J).GE.E0) GO TO 202
201   CONTINUE
202   CONTINUE
      IRM=IRANGE-1
      ED=E0-ZENGY(IRM)
C
C   REFLECTION PROBABILITY FOR FAST PARTICLE REFLECTION MODEL: RPROB
C
      IF (RINTG.GT.0.D0) THEN
        RPROB=MIN(PRFCT,RINTG)
      ELSE
        PRBRF=MIN(1.D0,MAX(0.D0,ZR(IRM)+QUOTR(IRM)*ED))
        RPROB=1.-(1.-PRBRF)*(COSIN**EXPP)
c slmod begin - tr
c        IF (PRFCF.LT.0.0) THEN
c          RPROB=MIN(-PRFCF,PRFCT)
c        ELSE
c          RPROB=MIN(RPROB*PRFCF,PRFCT)
c        ENDIF
c
        RPROB=MIN(RPROB*PRFCF,PRFCT)
c slmod end
      ENDIF
C
C  RELATIVE FRACTION OF COSINE VS. SPECULAR REFLECTION
      IF (EXPI.EQ.0.D0.OR.EXPI.GE.100.D0) THEN
        F1=1.
        F2=0.
      ELSE
        F1=MAX(0.015D0,MIN(1.D0,COSIN**EXPI))
        F2=SQRT(1.-F1*F1)
      ENDIF
C
      IATM=IGASF
      ISPZ=IATM
      WFAC=1.
      FR1=RANF_EIRENE( )
C
C  THERMAL PARTICLE MODEL
      IF (FR1.GE.RPROB) THEN
        IF (IGAST) 500,700,600
      ENDIF
C
C  FAST PARTICLE REFLECTION MODEL
C
C   ENERGY-REFLECTION COEFFICIENT
C   EPROB=1.-(1.-ERBRF)*(COSIN**EXPE)
C
      EFAC=COSIN**EXPE
      ESUM=E0-E0*EFAC
C
      IF (EINTG.GT.0.D0) THEN
        E0=E0*EINTG
C     ELSEIF (EINTG.LT.0.D0) THEN
C  E0 FROM MEAN ENERGY MODEL
C       E0=ESUM+(E0AV(IRM)+QUOTE(IRM)*ED)*EFAC
      ELSE
C  E0 FROM STOCHASTIC BEHRISCH MATRIX MODEL
C   REFLECTION ENERGY, "BEHRISCH MATRIX"
C   NUMBER OF BOXES IN THIS RANGE: IRANGE
C   DISTRIBUTION ZIDE(...,IRANGE)
C
        ZEP1=RANF_EIRENE( )
C
        DO 305 J=1,IRM
          IBOX=J
          ZDELTA=ZDE(J)
          ZE=ZRANGE(J)
          ZA=ZIDE(J,IRANGE)
          IF (ZA.GT.ZEP1) GO TO 307
305     CONTINUE
C  LAST BOX
        IBOX=IRANGE
        ZDELTA=ZDEL(IRANGE)
        ZE=ZENGY(IRANGE)
        ZA=1.
C
C   REFLECTION ENERGY, LINEAR INTERPOLATION
307     CONTINUE
        ZE0=ZE-(ZA-ZEP1)*ZDELTA/ZIDED(IBOX,IRANGE)
        ZE0=ZE0*E0/ZENGY(IRANGE)
        E0=ESUM+ZE0*EFAC
      ENDIF
C
C  E0 IS FOUND NOW. NEXT:
C  NEW WEIGHT, RESCALE ENERGY, SET VELOCITY
C
350   WEIGHT=WEIGHT*WFAC
      E0=E0/EREDUC
      VEL=RSQDVA(IATM)*SQRT(E0)
      ITYP=1
C     GOTO 400
C
400   CONTINUE
      IF (EXPI.LT.100.) THEN
        IF (F1.GT.0.999999) THEN
C  NO SPECULAR CONTRIBUTION (F2 = 0., F1 = 1.)
          IF (INIV4.LE.0) CALL FCOSIN
          VX=FC1(INIV4)
          VY=FC2(INIV4)
          VZ=FC3(INIV4)
          INIV4=INIV4-1
          CALL ROTATF (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ)
        ELSE
C  INCLUDE SPECULAR CONTRIBUTION (F2 > 0., F1 < 1.)
          ZTHET=PI2A*RANF_EIRENE( )
          ZSTHET=SIN(ZTHET)
          ZCTHET=COS(ZTHET)
          A=RANF_EIRENE( )
          ZCPHI=SQRT(A)
          ZSPHI=SQRT(1.-A)
C
          ZSPHI=ZSPHI*F1
          ZCPHI=SQRT(1.-ZSPHI*ZSPHI)
C
          VX=-ZCPHI*F1+        ZSPHI*ZSTHET*F2
          VY= ZSPHI*ZCTHET
          VZ= ZSPHI*ZSTHET*F1+ ZCPHI*F2
C
          CALL ROTATE (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ,COSIN)
        ENDIF
      ELSE
C   PURELY SPECULAR REFLECTION :EXPI .GE. 100 .
C   EXPI.GE.100 MEANS: INELASTIC+SPECULAR
        COSI2=-(COSIN+COSIN)
        VELX=VELX+COSI2*CRTX
        VELY=VELY+COSI2*CRTY
        VELZ=VELZ+COSI2*CRTZ
      ENDIF
      RETURN
C
C  "THERMAL MOLECULE MODEL"
C
C  CREATE MOLECULE OF SPECIES IMOL
C  WITH PROBABILITY PRFCT-RPROB, WHERE RPROB IS THE
C  PROBABILITY FOR BACKSCATTERING OF HOT ATOMS (.LE. PRFCT)
C  THE CONDITION FR1.GE.RPROB IS FULLFILLED AT THIS POINT
C
500   CONTINUE
      ITYP=2
      IMOL=-IGAST
      IF (IMOL.LT.1.OR.IMOL.GT.NMOLI) THEN
        FR2=RANF_EIRENE( )
        DO 501 I=1,NMOLI
          IMOL=I
          IF (FR2.LE.DMOL(IMOL)) GOTO 502
501     CONTINUE
        WRITE (6,*) 'SPECIES INDEX OUT OF RANGE IN REFLEC '
        WRITE (6,*) 'IMOL, MSURF ',IMOL,MSURF
        CALL EXIT
502     CONTINUE
      ENDIF
      ISPZ=NATMI+IMOL
C
C  FAST FRACTION
C     RPROBF=RPROB
C  THERMAL MOLECULE FRACTION
      RPROBM=PRFCT-RPROB
C  LOST FRACTION
      RPROBL=1.D0-PRFCT
C
      IF (WEIGHT.LT.WMIN.OR.RPROBM.LE.0.D0) THEN
C  NO SUPRESSION OF ABSORPTION
        PRTEST=RPROB+RPROBM
C  AT THIS POINT: 1.D0.GE.FR1.GE.RPROB
        IF (FR1.GT.PRTEST) GOTO 700
      ELSE
C  SUPRESSION OF ABSORPTION
        WMOLEC=RPROBM/(1.D0-RPROB)
        WLOSS =RPROBL/(1.D0-RPROB)
        IF (WLOSS.GT.0.D0) THEN
          WABS=WEIGHT*WLOSS
          IF (MSURF.GT.0) SPUMP(ISPZO,MSURF)=SPUMP(ISPZO,MSURF)+WABS
        ENDIF
        WEIGHT=WEIGHT*WMOLEC
      ENDIF
C
C  NUMBER OF MOLECULES PER INCIDENT PARTICLE
C  NOTE: ABSORPTION DUE TO RECOMBINATION OF ATOMS (ONLY A FRACTION OF
C  A MOLECULE IS RE-EMITTED PER INCIDENT ATOM) IS ALWAYS SUPRESSED
C
      FLPRT=DBLE(NPRIN)/DBLE(NPRT(ISPZ))
      WEIGHT=WEIGHT*FLPRT
C
C  REFLECT THERMAL MOLECULE
      IF (E0TERM.GT.0.D0) THEN
C  MONOENERGETIC, E0 (EV),  COSINE
        E0=E0TERM
        VEL=RSQDVM(IMOL)*SQRT(E0)
        F1=1.
        GOTO 400
      ELSEIF (E0TERM.LT.0.D0) THEN
C  SAMPLE FROM MAXWELLIAN FLUX AROUND INNER (!) NORMAL AT TEMP. TW (EV)
        TW=-E0TERM
        CALL VELOCS (TW,0.D0,0.D0,0.D0,0.D0,0.D0,RSQDVM(IMOL),
     .                CVRSSM(IMOL),
     .               -CRTX,-CRTY,-CRTZ,
     .               E0,VELX,VELY,VELZ,VEL)
      ELSE
        GOTO 991
      ENDIF
      RETURN
C
C  "THERMAL ATOM MODEL"
C
C    ONE ATOM IS BORN,
C    WITH PROBABILITY PRFCT-RPROB, WHERE RPROB IS THE
C    PROBABILITY FOR BACKSCATTERING OF HOT ATOMS (.LE. PRFCT)
C    THE CONDITION FR1.GE.RPROB IS FULLFILLED AT THIS POINT
C
C    E0TERM > 0: COSINE DISTRIBUTED, E0=E0TERM,
C    E0TERM < 0: MAXWELLIAN AT TEMP. T=-E0TERM
C    E0TERM = 0: THOMPSON DISTRIBUTION WITH SURF. BIND. EN.= EBIND
C                UNTIL NOW: ONLY FOR THERMAL ATOM REFLECTION MODEL
C
600   CONTINUE
      ITYP=1
      IATM=IGAST
      IF (IATM.GT.NATMI) THEN
        FR2=RANF_EIRENE( )
        DO 610 I=1,NATMI
          IATM=I
          IF (FR2.LE.DATM(IATM)) GOTO 611
610     CONTINUE
        WRITE (6,*) 'SPECIES INDEX OUT OF RANGE IN REFLEC '
        WRITE (6,*) 'IATM, MSURF ',IATM,MSURF
        CALL EXIT
611     CONTINUE
      ENDIF
      ISPZ=IATM
C
C  FAST FRACTION
C     RPROBF=RPROB
C  THERMAL ATOM FRACTION
      RPROBA=PRFCT-RPROB
C  LOST FRACTION
      RPROBL=1.D0-PRFCT
C
      IF (WEIGHT.LT.WMIN.OR.RPROBA.LE.0.D0) THEN
C  NO SUPRESSION OF ABSORPTION
        PRTEST=RPROB+RPROBA
C  AT THIS POINT: 1.D0.GE.FR1.GE.RPROB
        IF (FR1.GE.PRTEST) GOTO 700
      ELSE
C  SUPRESSION OF ABSORPTION
        WATOM=RPROBA/(1.D0-RPROB)
        WLOSS=RPROBL/(1.D0-RPROB)
        IF (WLOSS.GT.0.D0) THEN
          WABS=WEIGHT*WLOSS
          IF (MSURF.GT.0) SPUMP(ISPZO,MSURF)=SPUMP(ISPZO,MSURF)+WABS
        ENDIF
        WEIGHT=WEIGHT*WATOM
      ENDIF
C
C  REFLECT THERMAL ATOM
      IF (E0TERM.GT.0.D0) THEN
C  MONOENERGETIC, E0 (EV), +  STANDARD, COSINE LIKE
        E0=E0TERM
        E0_MEAN=E0TERM
        VEL=RSQDVA(IATM)*SQRT(E0)
        F1=1.
        GOTO 400
      ELSEIF (E0TERM.LT.0.D0) THEN
C  SAMPLE FROM MAXWELLIAN FLUX AROUND INNER (!) NORMAL AT TEMP. TW (EV)
        TW=-E0TERM
        CALL VELOCS (TW,0.D0,0.D0,0.D0,0.D0,0.D0,RSQDVA(IATM),
     .                CVRSSA(IATM),
     .               -CRTX,-CRTY,-CRTZ,
     .               E0,VELX,VELY,VELZ,VEL)
        RETURN
      ELSEIF (E0TERM.EQ.0.D0) THEN
C  SAMPLE FROM ENERGY FROM THOMPSON DISTRIBUTION + STAND. ANGULAR DISTR.
        E0=THOMP(EBIND,EMAXR)
        VEL=RSQDVA(IATM)*SQRT(E0)
        F1=1.
        GOTO 400
      ENDIF
C
C  ABSORB PARTICLE AT THIS SURFACE
C
700   CONTINUE
      IF (MSURF.GT.0) SPUMP(ISPZO,MSURF)=SPUMP(ISPZO,MSURF)+WEIGHT
      LGPART=.FALSE.
      WEIGHT=0.
      ITYP=0
      RETURN
C
C  ERROR MESSAGES FROM SUBR. REFLEC
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. REFLEC '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (6,*) 'MSURF = ',MSS
      WRITE (6,*) 'COSIN.LT.0. ', COSIN
      WRITE (6,*) 'STOP HISTORY NO. NPANU= ',NPANU
      GOTO 995
C
991   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. REFLEC '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (6,*) 'MSURF = ',MSS
      WRITE (6,*) 'E0TERM=0 '
      WRITE (6,*) 'STOP HISTORY NO. NPANU= ',NPANU
C
995   IF (NLTRC)  CALL CHCTRC(X0,Y0,Z0,16,15)
      LGPART=.FALSE.
      WEIGHT=0.
      RETURN
      END
