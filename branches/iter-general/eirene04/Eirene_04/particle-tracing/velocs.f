C
      SUBROUTINE VELOCS (TIWL,ESHET,VWL,VXWL,VYWL,VZWL,RSQDV,CVRSS,
     .                   CX,CY,CZ,
     .                   E0S,VELXS,VELYS,VELZS,VELS)
C
C  FETCH A NEW VELOCITY FROM A MAXWELLIAN FLUX AT A SURFACE GIVEN
C  BY THE NORMAL: CX,CY,CZ
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
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CRAND
      USE COMPRT
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: TIWL, ESHET, VWL, VXWL, VYWL, VZWL,
     .                      RSQDV, CVRSS, CX, CY, CZ
      REAL(DP), INTENT(OUT) :: E0S, VELXS, VELYS, VELZS, VELS
      REAL(DP) :: ARBV, A1, A2, A3, VLLX, VLLY, VLLZ, VMX, SHIFT,
     .          CCM, A4, FNOM, VMXSQ, FACTOR, VELSQ, VELSH, A5, A6,
     .          VFKT, VLX, RCCM, ZARG, ZARG2, VXDR, VYDR, VZDR, ERF
      REAL(DP), EXTERNAL :: RANF_EIRENE
C
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