C
      SUBROUTINE VELOCX(K,VXO,VYO,VZO,VLO,IOLD,NOLD,VELQ,NFLAG,
     .                  IRCX,DUMT,DUMV)
C
C  THIS SUBROUTINE CARRIES OUT AN ELASTIC COLLISION OF A TEST PARTICLE
C  WITH A BULK PARTICLE.
C  IT RETURNS THE POST COLLISION VELOCITY VECTOR.
C
C  NFLAG= 1:       SAMPLING FROM MONOENERGETIC DISTRIBUTION
C                  OF ION SPEED IN 1D, X DIRECTION
C                               IN 2D, X,Y DIRECTION
C                               IN 3D, X,Y,Z DIRECTION
C                  (I.E., DELTA FUNCTION IN ENERGY SPACE)
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
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CRAND
      USE CINIT
      USE CZT1
      USE CTRCEI
      USE COMPRT
      USE COMXS
      USE CLAST

      IMPLICIT NONE         

      REAL(DP), INTENT(IN) :: DUMT(3), DUMV(3)
      REAL(DP), INTENT(IN) :: VXO, VYO, VZO, VLO
      REAL(DP), INTENT(OUT) :: VELQ
      INTEGER, INTENT(IN) :: K, IOLD, NOLD, NFLAG, IRCX

      REAL(DP) :: VXN, VYN, VZN, VX,VY,VZ, VN, ZARGX, ZARGY, ZARGZ,
     .          VXDR, VYDR, VZDR, VRELQ, E0MAX, TIMAX, SIGS, VRELS,
     .          ELABS, WRMEAN, TEST, VREL, WRAT, WO, ELAB, CXS, 
     .          VR, VRQ, CROSS, ELMAX, ELMIN
      REAL(DP), EXTERNAL :: RANF_EIRENE
      REAL(DP) :: WRMAX = -1000.D0, WRMIN = 1000.D0

      INTEGER :: IMEAN, IRX, ICOUNT, J, JJ, IRL, IREAC
      INTEGER :: IPMAX = 0, IPMIN = 0, IFIRST = 0

      SAVE
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        DO IRL=1,NRCXI
          IFLRCX(IRL)=0
          NCMEAN(IRL)=0
          XCMEAN(IRL)=0.D0
        ENDDO
      ENDIF
C
      IF (IFLRCX(IRCX).EQ.0.AND.NFLAG.NE.2) THEN
        IFLRCX(IRCX)=-1
C  PREPARE REJECTION SAMPLING OF INCIDENT ION VELOCITY
C  IS CROSS SECTION AVAILABLE?
        IREAC=MODCOL(3,1,NOLD,IPLS)
        IF (IREAC.EQ.0) GOTO 1
C
        elmin=log(0.1)
        elmax=log(1.e4)
        SGCVMX(IRCX)=-1.D60
        JJ=1
        do j=1,1000
          elab=elmin+(j-1)/999.*(elmax-elmin)
          CXS=CROSS(ELAB,IREAC,IRCX,'VELOCX 1')
          vrq=exp(elab-defCX(IRCX))
          vr=sqrt(vrq)
          if (cXS*vr.gt.SGCVMX(IRCX)) then
            JJ=J
            SGCVMX(IRCX)=cXS*vr
          endif
        enddo
        CALL LEER(1)
        WRITE (6,*) 'FIRST CALL TO VELOCX FOR IRCX= ',IRCX
        WRITE (6,*) 'PREPARE REJECTION TECHNIQUE '
        WRITE (6,*) 'SGCVMX IN VELOCX,JJ ',SGCVMX(IRCX),JJ
        IF (JJ.NE.1.AND.JJ.NE.1000) IFLRCX(IRCX)=1
        CALL LEER(1)
      ENDIF
1     CONTINUE
C
      ICOUNT=1

      IF (K.GT.0) THEN
        ZARGX=ZRG(IPLS,K)
        ZARGY=ZRG(IPLS,K)
        ZARGZ=ZRG(IPLS,K)
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL VECUSR(2,VXDR,VYDR,VZDR,IPLS)
          ELSE
            VXDR=VXIN(IPLS,K)
            VYDR=VYIN(IPLS,K)
            VZDR=VZIN(IPLS,K)
          END IF
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
123   CONTINUE
      IF (INIV2.LE.0) CALL FGAUSS
C
C  SAMPLE FROM 3D MAXWELLIAN
      VXN=FG1(INIV2)
      VYN=FG2(INIV2)
      VZN=FG3(INIV2)
      INIV2=INIV2-1
C
C  DRIFTING, MONOENERGETIC ISOTROPIC DISTRIBUTION
C
      IF (NFLAG.EQ.1) THEN
C  ZT1 CORRESPONDS TO ROOT MEAN SQUARE VELOCITY AT TIIN(IPLS,K)
        VEL=SQRT(ZT1(IPLS,K))
        VN=VEL/SQRT(VXN*VXN+VYN*VYN+VZN*VZN)
        VXN=VXN*VN+VXDR
        VYN=VYN*VN+VYDR
        VZN=VZN*VN+VZDR
      ELSE
        VXN=VXN*ZARGX+VXDR
        VYN=VYN*ZARGY+VYDR
        VZN=VZN*ZARGZ+VZDR
      ENDIF
C
C  DRIFTING MAXWELLIAN DISTRIBUTION (FOR MAXWELL-POTENTIAL: SIGMA*V = CONST.)
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
        RETURN
C
      ELSE
C
        VX=VXO*VLO
        VY=VYO*VLO
        VZ=VZO*VLO
C
C   ALL OTHER DISTRIBUTIONS
C
C   WEIGHT CORRECTION DUE TO ENERGY DEPENDENCE IN CROSS SECTION
C   OR: REJECTION     DUE TO ENERGY DEPENDENCE IN CROSS SECTION
C   PRESENT VERSION: REJECTION
        VRELQ=(VXN-VX)**2+(VYN-VY)**2+(VZN-VZ)**2
        VREL=SQRT(VRELQ)
        ELAB=LOG(VRELQ)+DEFCX(IRCX)
        IREAC=MODCOL(3,1,NOLD,IPLS)
        CXS=CROSS(ELAB,IREAC,IRCX,'VELOCX 2')
C
C       IF (NLREJC) THEN
        IF (IFLRCX(IRCX).GT.0) THEN
          TEST=RANF_EIRENE()*SGCVMX(IRCX)
          if (test.gt.cxs*vrel) then
c  reject
            icount=icount+1
            if (icount.lt.500) goto 123
            write (6,*) 'icount too large IN VELOCX. ACCEPT SAMPLE '
            write (6,*) 'npanu, ireac, ircx, ELAB ',
     .                   npanu, ireac, ircx, ELAB
          else
c  accept
            xcmean(ircx)=xcmean(ircx)+icount
            ncmean(ircx)=ncmean(ircx)+1
          endif
C       ELSEIF (NLWEIGHT) THEN
        ELSE
          WEIGHT=WEIGHT*CXS*VREL*DIIN(IPLS,K)/SIGVCX(IRCX)
        ENDIF
C
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
      ENDIF
C
      RETURN
C
999   CONTINUE
      WRITE (6,*) 'PARAMETER ERROR IN SUBR. VELOCX. EXIT CALLED'
      CALL EXIT_OWN(1)
      END
