C  27.6.05: iadd removed
C  02.01.06:    SIGMAX NOW SET ONLY FOR ACTIVE REACTIONS
C               added: jcou,ncou
!pb  12.10.06:  modcol revised
!pb  28.11.06:  initialization of XSTOR reactivated because of trouble in
!pb             BGK iteration
C
      FUNCTION EIRENE_FPATHPH (K,CFLAG,JCOU,NCOU)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR PHOTON
C   "BEAM" OF VELOCITY (E0,VEL_X,Y,Z) IN DRIFTING MAXWELLIAN BACKGROUND
C   IN CELL K
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CZT1
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_PHOTON
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(OUT) :: CFLAG(7,3)
      INTEGER, INTENT(IN) :: K, JCOU,NCOU
 
      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLSV)
      REAL(DP) :: EIRENE_FPATHPH, sigmax, sigv, eirene_feplot3, 
     .            vx, vy, vz,
     .            DENEL, PVELQ0, fac
      integer :: il, kk, ipl, irot, ipot, j, iph, i1, i2
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTORV=0.D0
!pb      IF (NCOU.GT.1) THEN
        XSTOR=0.D0
!pb      ENDIF
      EIRENE_FPATHPH = 1.E10_DP
      SIGMAX=0.D0
 
      IF (LGVAC(K,0)) RETURN
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)
      DO 2 IPLS=1,NPLSI
        ZTI(IPLS)=ZT1(IPLS,K)
2       DENIO(IPLS)=DIIN(IPLS,K)
C
      PVELQ0=VEL*VEL
      DO 3 IPLS=1,NPLSV
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL EIRENE_VECUSR (2,VX,VY,VZ,IPLS)
          ELSE
            VX=VXIN(IPLS,K)
            VY=VYIN(IPLS,K)
            VZ=VZIN(IPLS,K)
          END IF
          PVELQ(IPLS)=(VELX*VEL-VX)**2+
     .                (VELY*VEL-VY)**2+
     .                (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLS)=PVELQ0
        ENDIF
3     CONTINUE
C
csw
csw OT processes (photonic reactions)
csw
60    CONTINUE
      if(phv_lgphot(iphot,0,0) == 0) goto 70
      do 61 ipot=1,phv_nphoti(iphot)
        irot=phv_lgphot(iphot,ipot,0)
        ipls =phv_lgphot(iphot,ipot,1)
        il   =phv_lgphot(iphot,ipot,2)
        kk   =phv_lgphot(iphot,ipot,3)
        IF (LGVAC(K,IPLS)) GOTO 61
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(7,2,   IROT).EQ.1) THEN
          GOTO 997
        ELSEIF (MODCOL(7,2,   IROT).EQ.2) THEN
C  MODEL 2:
C  BEAM - MAXWELLIAN RATE. FULL ACCOUNT FOR DOPPLER SHIFT
cdr       kk   = nreaot(irot)
cdr   effective energy e0_eff due to doppler shift from directed motion
cdr       e0_eff=
cdr  getcoeff liefert nun maxw. average ueber Ti(ipls), z.b. voigt, ....
          call EIRENE_PH_GETCOEFF(kk,iphot,0,k,ipls,fac,sigv)
          sigv=sigv*diin(ipls,k)
          if(phv_muldens .EQ. 0) then
cdr  hier in fpathph kann es keine spontanen raten (1/s) geben.
cdr  spaeter: allgemein raten (1/s) auch fuer testteilchen (fpatha, fpathm, fpat
cdr           als neue option einfuehren, analog Aik in xsectp.
            GOTO 997
          endif
          SIGVOT(irot)=sigv
          GOTO 997
        ELSEIF (MODCOL(7,2,   IROT).EQ.4) THEN
C  MODEL 4:
C  BEAM-BEAM RATE. IGNORE DOPPLER SHIFT DUE TO THERMAL MOTION,
C                  INCLUDE DOPPLER SHIFT DUE TO DIRECTED MOTION
cdr       kk   = nreaot(irot)
cdr   effective energy e0_eff due to doppler shift from directed motion
cdr       e0_eff=
cdr       ireac=modcol(7,1,irot)
cdr  ireac entspricht "typ" in Getcoeff - cross section (lorentz, vdw, ...)
cdr  allerdings kann hier der "querschnitt" von hintergrundparametern abhaengen
          call EIRENE_PH_GETCOEFF(kk,iphot,0,k,ipls,fac,sigv)
          sigv=sigv*diin(ipls,k)
cdr  besser: sig vc. E0_effective, ohne v = vrel = c, dann
cdr  dann:   sigv=sig * c * diin
cdr  denn:   sigv enthaelt hier keine faltung ueber background maxw. vel-verteil
          if(phv_muldens .EQ. 0) then
cdr  hier in fpathph kann es keine spontanen raten (1/s) geben.
cdr  spaeter: allgemein raten (1/s) auch fuer testteilchen (fpatha, fpathm, fpat
cdr           als neue option einfuehren, analog Aik in xsectp.
            GOTO 997
          endif
          SIGVOT(irot)=sigv
cdr       ESIGOT(irot,1)=e0*sigv   ziemlich sicher falsch
        ELSE
          GOTO 997
        ENDIF
 
        SIGMAX=MAX(SIGMAX,SIGVOT(IROT))
        SIGOTT=SIGOTT+SIGVOT(IROT)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(7,4,     IROT).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGOT(IROT,1)=EPLOT3(IROT,K,1)
          ELSE
            ESIGOT(IROT,1)=EIRENE_FEPLOT3(IROT,K)
          END IF
          CFLAG(7,1)=2
        ELSEIF (MODCOL(7,4,     IROT).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGOT(IROT,1)=EPLOT3(IROT,K,1)
          ELSE
            ESIGOT(IROT,1)=EIRENE_FEPLOT3(IROT,K)
          END IF
          CFLAG(7,1)=1
        ELSE
          GOTO 997
        ENDIF
61    CONTINUE
 
70    CONTINUE
c
C     TOTAL
C
100   CONTINUE
 
C
      IF (SIGOTT.GT.0._DP) THEN
        DO IROT=1,NROT
          IF (SIGVOT(IROT) .LE. SIGMAX*1.D-10) THEN
            SIGOTT=SIGOTT-SIGVOT(IROT)
            SIGVOT(IROT) = 0.D0
          END IF
        END DO
      END IF
 
      SIGTOT=SIGEIT+SIGPIT+SIGCXT+SIGELT+SIGOTT
      IF (SIGTOT.GT.1.D-20) THEN
        EIRENE_FPATHPH=VEL/SIGTOT
        ZMFPI=1./EIRENE_FPATHPH
      ENDIF
C
      RETURN
997   CONTINUE
      WRITE (iunout,*)
     .  'ERROR IN FPATHPH: INCONSISTENT PHOTON COLL. DATA.'
      WRITE (iunout,*) 'IPHOT,IROT,MODCOL(7,J,IROT) '
      WRITE (iunout,*) IPHOT,IROT,(MODCOL(7,J,IROT),J=1,4)
      CALL EIRENE_EXIT_OWN(1)
      END
