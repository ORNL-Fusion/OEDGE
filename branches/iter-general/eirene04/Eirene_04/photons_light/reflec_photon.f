C
C
      SUBROUTINE REFLEC_photon
C
C  REFLECT ESCAPING PHOTONS
C  INPUT:
C       ILREF = 0  PERFECT ABSORPTION
C       ILREF = 1  ERIC HOLLMAN DATABASE
C
C       ITYP  = 0  INCIDENT PHOTON
C  OUTPUT:
C     LGPART= TRUE AND:
C       ITYP = 0  PHOTON IPHOT IS RETURNED TO CALLING PROGRAM
C     LGPART= FALSE  NO PARTICLE IS RETURNED (ABSORBTION)
C       ITYP = -1
C
      USE PRECISION
      USE PARMMOD
      USE PHOT_REFLEC
      USE CTRCEI
      USE CADGEO
      USE CLGIN
      USE COMUSR
      USE CLOGAU
      USE CESTIM
      USE COMPRT
      USE CRAND
      USE CREF
      USE CCONA
      USE PHOTON
      USE CPES

      IMPLICIT NONE
    
      REAL(DP), INTENT(IN) :: WMIN, XMP, XCP
      INTEGER, INTENT(IN) :: NPRIN
      INTEGER, INTENT(INOUT) :: IGASF, IGAST
      INTEGER, SAVE :: IFIRST=0, NPANOLD=0
      INTEGER :: ILIM, ISP, ICOUNT, ISTS, ISEE, MODREF, ISPZO,
     .           IMAT, IREFL, MSS
      REAL(DP) :: DUMMY, XMW, XCW, E0TERM, EBIND, PRFCF, PRFCT,
     .            EXPP, EXPE, EXPI, RINTG, EINTG, AINTG, COSIN,
     .            THETA, XLAMBDA, THETA_OUT, ALPHA_OUT, WFAC, 
     .            FR1, ZCPHI, ZSPHI, ZCTHET, ZSTHET, VX, VY, VZ,
     .            RPROB,WABS
      REAL(DP), EXTERNAL :: RANF_EIRENE, ranset_eirene
      INTEGER, EXTERNAL :: RANGET_EIRENE, IDEZ
C
C---------------------------------------------------------------------
C

C
C  INITIALIZE SURFACE REFLECTION MODELS FOR PHOTONS
C
      ENTRY REFLC0_PHOTON
C
      IF (IFIRST.EQ.1) RETURN
      IFIRST=1
C
      CALL INIT_REFL_HLM
C
      IF (TRCREF) THEN
        CALL LEER(2)
      ENDIF
C
C
C  PRINTOUT REFLECTION PROPERTIES OF SURFACES
C  ALREADY DONE FROM REFLC0
C
C
      RETURN
C
      ENTRY REFLC1_PHOTON (WMIN,XMP,XCP,NPRIN,IGASF,IGAST)
C
C  SYNCHRONIZE RANDOM NUMBERS
C
      IF (NLCRR.AND.(NPANU.NE.NPANOLD)) THEN
C  INITIALIZE RANDOM NUMBERS FOR EACH PARTICLE, TO GENERATE CORRELATION
C       Call RANSET_EIRENE(ISEED)
        dummy=ranset_eirene(iseedR)
        DUMMY=RANF_EIRENE( )
        ISEEDR=ranget_eirene(isee)
        ISEEDR=INTMAX-ISEEDR
        NPANOLD=NPANU
      END IF
C
C  SURFACE NUMBER  : MSURF (MSURF=0: DEFAULT MODEL)
C  SPECIES INDEX   : ISPZ
C
      MODREF=IDEZ(ILREF(MSURF),2,2)
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
C
C
C   TENTATIVELY ASSUME  REFLECTION
      LGPART=.TRUE.
C   COSINE OF ANGLE OF INCIDENCE
      COSIN=VELX*CRTX+VELY*CRTY+VELZ*CRTZ
      IF (COSIN.LT.0.D0) GOTO 993
C
C
C
C   MODREF=0: "PERFECTLY ABSORBING SURFACE, DEFAULT
C   MODREF=1: "DATABASE REFLECTION MODEL" 
C
      IF (MODREF.EQ.1.AND.PRFCF.GT.0.) THEN
        GOTO 100
      ELSE
C  ABSORB THIS PHOTON
        GOTO 700
      ENDIF
C
C  DATABASE REFLECTION MODEL STARTS HERE
C
100   CONTINUE
C
C   CHECK IF WALL REFLECTION DATA FOR IPHOT INCIDENT ON
C   XWALL/ZWALL ARE AVAILABLE
C
C  CARBON OR MOLYBDENUM ?

      IMAT = 2
      
C
      THETA = ACOS(COSIN)*RADDEG
      XLAMBDA = hpcl/E0*10._DP*1.E7_DP

      CALL REFLECT_HOLLMANN (THETA, XLAMBDA, IMAT, IREFL, THETA_OUT,
     .                       ALPHA_OUT, RPROB)
C
130   CONTINUE
C
C
C   DECIDE IF PARTICLE IS TO BE REFLECTED OR ABSORBED
C   (NO THERMAL RE-EMISSION MODEL FOR INCIDENT PHOTONS)
C
      IF (WEIGHT.GT.WMIN) THEN
C  WITH SUPPRESSION OF ABSORPTION
        WABS=WEIGHT*(1.D0-RPROB)
        IF (WABS.GT.0.D0) THEN
          IF (LSPUMP) SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WABS
        ENDIF  
        WEIGHT=WEIGHT-WABS
        IF (WEIGHT.LE.EPS30) GOTO 700
      ELSE
C  NO SUPPRESSION OF ABSORPTION
        FR1=RANF_EIRENE( )
        IF (FR1.GT.RPROB) GOTO 700
      ENDIF
C
C  SPECIES OF REFLECTED PARTICLE
      IF (IGASF.LT.1.OR.IGASF.GT.NPHOTI) GOTO 992
      IPHOT=IGASF
      ISPZ=IPHOT
      ITYP=0
C
C  ENERGIE (WAVELENGTH):  NOT MODIFIED
C
C
C  POLAR ANGLE OF REFLECTION  (THETA BEI OLIVER)
C
C
      ZCPHI=COS(THETA_OUT)
C  LIMIT COSINE OF POLAR ANGLE TO 85. DEGREES
C  (I.E., 5 DEGREES AGAINST SURFACE TANGENTIAL PLANE)
      ZCPHI=MIN(0.999999_DP,MAX(0.08716_DP,ZCPHI))
      ZSPHI=SQRT(1.-ZCPHI*ZCPHI)
C
C  AZIMUTAL ANGLE OF REFLECTION    (ALPHA BEI OLIVER)
C
      ZCTHET=COS(ALPHA_OUT)
      ZCTHET=MAX(-.999999_DP,MIN(0.999999_DP,ZCTHET))
      ZSTHET=SQRT(1.-ZCTHET*ZCTHET)
      ZSTHET=ZSTHET*SIGN(1._DP,(RANF_EIRENE( )-0.5_DP))
C
      VX=-ZCPHI
      VY=ZSPHI*ZSTHET
      VZ=ZSPHI*ZCTHET
      IF (COSIN.GT.0.999999) THEN
        CALL ROTATF (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ)
      ELSE
        CALL ROTATE (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ,COSIN)
      ENDIF
      RETURN
C
C  ABSORB PARTICLE AT THIS SURFACE
C
700   CONTINUE
      IF (LSPUMP.AND.(MSURF.GT.0)) 
     .   SPUMP(ISPZO,MSURF)=SPUMP(ISPZO,MSURF)+WEIGHT
      LGPART=.FALSE.
      WEIGHT=0.
      ITYP=-1
      RETURN
C
C  ERROR MESSAGES FROM SUBR. REFLEC_PHOTON
C
C
992   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. REFLEC_PHOTON '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (6,*) 'MSURF = ',MSS
      WRITE (6,*) 'IGASF, IGAST ?? '
      WRITE (6,*) 'STOP HISTORY NO. NPANU= ',NPANU
      GOTO 999
c
993   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. REFLEC_PHOTON '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (6,*) 'MSURF = ',MSS
      WRITE (6,*) 'COSIN.LT.0. ', COSIN
      WRITE (6,*) 'STOP HISTORY NO. NPANU= ',NPANU
      GOTO 999
C
999   IF (NLTRC)  CALL CHCTRC(X0,Y0,Z0,16,18)
      LGPART=.FALSE.
      WEIGHT=0.
      RETURN
C
      END