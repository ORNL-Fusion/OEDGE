C
      SUBROUTINE VELOEI(K,IRDS,VXO,VYO,VZO,VLO)
C
C  FETCH A NEW VELOCITY OF TEST PARTICLE AFTER ELECTRON IMPACT COLLISION
C
C  K   : CELL INDEX
C  VXO : X COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VYO : Y COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VZO : Z COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VLO : VELOCITY OF TEST PARTICLE BEFORE EVENT
C
C
C  FIND TYPE OF NEXT GENERATION PARTICLE (ATOM, MOLECULE, TEST ION)
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CRAND
      USE CZT1
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: VXO, VYO, VZO, VLO
      INTEGER, INTENT(IN) :: K, IRDS
      REAL(DP) :: FEHVDS1, VXDIS, VYDIS, VZDIS, EHEAVY, VX, VY, VZ, 
     .            VELQ, CVRSS, RSQDV, EDISS, ZEP3, VELDS
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: ISPZI, ISPZM, ISPZA

      ZEP3=RANF_EIRENE( )
C
      IF (ZEP3.LE.P2ND(IRDS,NSPA)) THEN
C
C  A NEUTRAL ATOM IS BORN, FIND SPECIES INDEX IATM AND WEIGHT
C
        ITYP=1
        DO 448 IATM=1,NATMIM
          ISPZA=NSPH+IATM
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
          ISPZM=NSPA+IMOL
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
        CALL EXIT_OWN(1)
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
