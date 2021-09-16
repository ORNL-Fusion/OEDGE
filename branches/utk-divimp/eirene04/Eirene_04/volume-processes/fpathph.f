C
C
      FUNCTION FPATHPH (K,CFLAG)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR PHOTON
C   "BEAM" OF VELOCITY (E0,VEL_X,Y,Z) IN DRIFTING MAXWELLIAN BACKGROUND
C   IN CELL K
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CZT1
      USE COMPRT
      USE COMXS
      USE CSPEI
      USE PHOTON

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(7,3)
      INTEGER, INTENT(IN) :: K

      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLSV)
      REAL(DP) :: FPATHPH, sigmax, sigv, feplot3, vx, vy, vz,
     .            DENEL, PVELQ0, fac
      integer :: il, kk, ipl, irot, ipot, iadd, j, iph, i1, i2
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTOR=0.D0
      XSTORV=0.D0
c     DO I2=1,MSTOR2
c       DO I1=1,MSTOR1
c         XSTOR(I1,I2) = 0._DP
c       END DO
c     END DO
c     DO I1=1,NSTORV
c       XSTORV(I1) = 0._DP
c     END DO
      FPATHPH = 1.E10_DP
      iadd = phv_nrota

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
            CALL VECUSR (2,VX,VY,VZ,IPLS)
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
        IF (MODCOL(7,2,   IPHOT,IPLS).EQ.1) THEN
          GOTO 997
        ELSEIF (MODCOL(7,2,   IPHOT,IPLS).EQ.2) THEN
C  MODEL 2:
C  BEAM - MAXWELLIAN RATE. FULL ACCOUNT FOR DOPPLER SHIFT
cdr       kk   = nreaot(irot)
cdr   effective energy e0_eff due to doppler shift from directed motion
cdr       e0_eff=
cdr  getcoeff liefert nun maxw. average ueber Ti(ipls), z.b. voigt, ....
!pb          sigv = PH_GETCOEFF(kk,iphot,0,k,ipls,irot)
          call PH_GETCOEFF(kk,iphot,0,k,ipls,irot,fac,sigv)
          sigv=sigv*diin(ipls,k)
          if(phv_muldens .EQ. 0) then
cdr  hier in fpathph kann es keine spontanen raten (1/s) geben.
cdr  spaeter: allgemein raten (1/s) auch fuer testteilchen (fpatha, fpathm, fpathi...)
cdr           als neue option einfuehren, analog Aik in xsectp.
            GOTO 997
          endif
          SIGVOT(iadd+irot)=sigv
          GOTO 997
        ELSEIF (MODCOL(7,2,   IPHOT,IPLS).EQ.4) THEN
C  MODEL 4:
C  BEAM-BEAM RATE. IGNORE DOPPLER SHIFT DUE TO THERMAL MOTION,
C                  INCLUDE DOPPLER SHIFT DUE TO DIRECTED MOTION
cdr       kk   = nreaot(irot)
cdr   effective energy e0_eff due to doppler shift from directed motion
cdr       e0_eff=
cdr       ireac=modcol(7,1,iphot,ipls)
cdr  ireac entspricht "typ" in Getcoeff - cross section (lorentz, vdw, ...)
cdr  allerdings kann hier der "querschnitt" von hintergrundparametern abhaengen
!pb          sigv = PH_GETCOEFF(kk,iphot,0,k,ipls,irot)
          call PH_GETCOEFF(kk,iphot,0,k,ipls,irot,fac,sigv)
          sigv=sigv*diin(ipls,k)
cdr  besser: sig vc. E0_effective, ohne v = vrel = c, dann
cdr  dann:   sigv=sig * c * diin
cdr  denn:   sigv enthaelt hier keine faltung ueber background maxw. vel-verteilung
          if(phv_muldens .EQ. 0) then
cdr  hier in fpathph kann es keine spontanen raten (1/s) geben.
cdr  spaeter: allgemein raten (1/s) auch fuer testteilchen (fpatha, fpathm, fpathi...)
cdr           als neue option einfuehren, analog Aik in xsectp.
            GOTO 997
          endif
          SIGVOT(iadd+irot)=sigv
cdr       ESIGOT(iadd+irot,1)=e0*sigv   ziemlich sicher falsch
        ELSE
          GOTO 997
        ENDIF
        SIGOTT=SIGOTT+SIGVOT(iadd+irot)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(7,4,     IPHOT,IPLS).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGOT(IADD+IROT,1)=EPLOT3(IADD+IROT,K,1)
          ELSE
            ESIGOT(IADD+IROT,1)=FEPLOT3(IADD+IROT,K)
          END IF
          CFLAG(7,1)=2
        ELSEIF (MODCOL(7,4,     IPHOT,IPLS).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGOT(IADD+IROT,1)=EPLOT3(IADD+IROT,K,1)
          ELSE
            ESIGOT(IADD+IROT,1)=FEPLOT3(IADD+IROT,K)
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
csw extend 1:4 to 1:5 (OT process)
!      SIGMAX=MAXVAL(XSTOR(:,1:5))
!      WHERE (XSTOR(:,1:5) .LE. SIGMAX*1.D-10 )
!        XSTOR(:,1:5) = 0.D0
!      END WHERE
      SIGMAX=MAXVAL(XSTOR(1:mstor1,1:5))
      WHERE (XSTOR(1:mstor1,1:5) .LE. SIGMAX*1.D-10 )
        XSTOR(1:mstor1,1:5) = 0.D0
      END WHERE
c     DO I1=1,MSTOR1
c       IF (XSTOR(I1,1) .LE. SIGMAX*1.D-10 ) XSTOR(I1,1) = 0._DP
c       IF (XSTOR(I1,2) .LE. SIGMAX*1.D-10 ) XSTOR(I1,2) = 0._DP
c       IF (XSTOR(I1,3) .LE. SIGMAX*1.D-10 ) XSTOR(I1,3) = 0._DP
c       IF (XSTOR(I1,4) .LE. SIGMAX*1.D-10 ) XSTOR(I1,4) = 0._DP
c       IF (XSTOR(I1,5) .LE. SIGMAX*1.D-10 ) XSTOR(I1,5) = 0._DP
c     END DO
C
      SIGTOT=SIGEIT+SIGPIT+SIGCXT+SIGELT+SIGOTT
      IF (SIGTOT.GT.1.D-20) THEN
        FPATHPH=VEL/SIGTOT
        ZMFPI=1./FPATHPH
      ENDIF
C     
      RETURN
990   CONTINUE
      IPH=     IPHOT
      WRITE (6,*) 'ERROR IN FPATHPH: INCONSISTENT ELEC. IMP. DATA'
      WRITE (6,*) 'IPH,MODCOL(1,J,IPH,1) '
      WRITE (6,*) IPHOT,(MODCOL(1,J,IPH,1),J=1,4)
      CALL EXIT_OWN(1)
991   CONTINUE
      IPH=     IPHOT
      WRITE (6,*) 'ERROR IN FPATHPH: INCONSISTENT ION IMP. DATA'
      WRITE (6,*) 'IPH,IPL,MODCOL(4,J,IPH,IPL) '
      DO 993 IPL=1,NPLSI
        WRITE (6,*) IPHOT,IPL,(MODCOL(4,J,IPH,IPL),J=1,4)
993   CONTINUE
      CALL EXIT_OWN(1)
992   CONTINUE
      IPH=IPHOT
      WRITE (6,*) 'ERROR IN FPATHPH: INCONSISTENT CHARGE EXCHANGE DATA'
      WRITE (6,*) 'IPH,IPL,MODCOL(3,J,IPH,IPL) '
      DO 994 IPL=1,NPLSI
        WRITE (6,*) IPHOT,IPL,(MODCOL(3,J,IPH,IPL),J=1,4)
994   CONTINUE
      CALL EXIT_OWN(1)
995   CONTINUE
      IPH=IPHOT
      WRITE (6,*) 'ERROR IN FPATHPH: INCONSISTENT ELASTIC COLL. DATA'
      WRITE (6,*) 'IPH,IPL,MODCOL(5,J,IPH,IPL) '
      DO 996 IPL=1,NPLSI
        WRITE (6,*) IPHOT,IPL,(MODCOL(5,J,IPH,IPL),J=1,4)
996   CONTINUE
      CALL EXIT_OWN(1)
997   CONTINUE
      IPH=IPHOT
      WRITE (6,*) 'ERROR IN FPATHPH: INCONSISTENT PHOTON COLL. DATA'
      WRITE (6,*) 'IPH,IPL,MODCOL(7,J,IPH,IPL) '
      DO 998 IPL=1,NPLSI
        WRITE (6,*) IPHOT,IPL,(MODCOL(7,J,IPH,IPL),J=1,4)
998   CONTINUE
      CALL EXIT_OWN(1)
      END
