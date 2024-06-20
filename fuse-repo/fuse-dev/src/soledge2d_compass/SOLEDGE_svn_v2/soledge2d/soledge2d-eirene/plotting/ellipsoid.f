 
 
      SUBROUTINE EIRENE_ELLIPSOID
     .  (X0,Y0,Z0,CX,CY,CZ,XLIMS1,YLIMS1,ZLIMS1,
     .                      XLIMS2,YLIMS2,ZLIMS2,RLB,ILCOL,NX,NY,NZ)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_CCONA
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: X0, Y0, Z0, CX, CY, CZ, XLIMS1, YLIMS1,
     .                      ZLIMS1, XLIMS2, YLIMS2, ZLIMS2, RLB
      INTEGER, INTENT(IN) :: ILCOL, NX, NY, NZ
      REAL(DP) :: XP(2*101), YP(2*101), PHIAN(9), PHIEN(9)
      REAL(DP) :: RX, RY, RZ, DPH, DZ, X, Y, Z, ZE, RAD, XS, YS, PHI,
     .          DX, DY, XA, XE, YA, YE, ZA
      INTEGER :: IPART, IPRT, N, IZ, IX, IY, I
      INTEGER :: IP(2*101)
 
      XA = MAX(X0-CX,XLIMS1)
      XE = MIN(X0+CX,XLIMS2)
      IF (NX > 1) THEN
        DX = (XE-XA)/FLOAT(NX-1)
      ELSE
        DX = 0.D0
      END IF
 
      YA = MAX(Y0-CY,YLIMS1)
      YE = MIN(Y0+CY,YLIMS2)
      IF (NY > 1) THEN
        DY = (YE-YA)/FLOAT(NY-1)
      ELSE
        DY = 0.D0
      END IF
 
      ZA = MAX(Z0-CZ,ZLIMS1)
      ZE = MIN(Z0+CZ,ZLIMS2)
      IF (NZ > 1) THEN
        DZ = (ZE-ZA)/FLOAT(NZ-1)
      ELSE
        DZ = 0.D0
      END IF
 
      N = 101
 
!  PLOT Z-ISOLINES
      CALL GRNWPN(ILCOL)
 
      DO IZ = 1, NZ
        Z = ZA + (IZ-1)*DZ
        RAD = 1.D0 - (Z-Z0)**2/CZ**2
        IF (RAD < 0.D0) CYCLE
        RX = CX * SQRT(RAD)
        RY = CY * SQRT(RAD)
 
        CALL EIRENE_CTELL (X0,Y0,RX,RY,XLIMS1,YLIMS1,XLIMS2,YLIMS2,RLB,
     .              PHIAN,PHIEN,IPART)
 
        DO IPRT=1,IPART
          DPH = (PHIEN(IPRT)-PHIAN(IPRT)) / FLOAT(N-1)
          X = X0 + RX*COS(PHIAN(IPRT))
          Y = Y0 + RY*SIN(PHIAN(IPRT))
          CALL EIRENE_PL3D (X,Y,Z,XS,YS)
          CALL GRJMP (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          DO I = 2,N
            PHI = PHIAN(IPRT) + (I-1)*DPH
            X = X0 + RX*COS(PHI)
            Y = Y0 + RY*SIN(PHI)
            CALL EIRENE_PL3D (X,Y,Z,XS,YS)
            CALL GRDRW (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          END DO
        END DO
 
      END DO ! IZ
 
!  PLOT Y-ISOLINES
 
      DO IY = 1, NY
        Y = YA + (IY-1)*DY
        RAD = 1.D0 - (Y-Y0)**2/CY**2
        IF (RAD < 0.D0) CYCLE
        RX = CX * SQRT(RAD)
        RZ = CZ * SQRT(RAD)
 
        CALL EIRENE_CTELL (X0,Z0,RX,RZ,XLIMS1,ZLIMS1,XLIMS2,ZLIMS2,RLB,
     .              PHIAN,PHIEN,IPART)
 
        DO IPRT=1,IPART
          DPH = (PHIEN(IPRT)-PHIAN(IPRT)) / FLOAT(N-1)
          X = X0 + RX*COS(PHIAN(IPRT))
          Z = Z0 + RZ*SIN(PHIAN(IPRT))
          CALL EIRENE_PL3D (X,Y,Z,XS,YS)
          CALL GRJMP (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          DO I = 2,N
            PHI = PHIAN(IPRT) + (I-1)*DPH
            X = X0 + RX*COS(PHI)
            Z = Z0 + RZ*SIN(PHI)
            CALL EIRENE_PL3D (X,Y,Z,XS,YS)
            CALL GRDRW (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          END DO
        END DO
      END DO ! IY
 
!  PLOT X-ISOLINES
 
      DO IX = 1, NX
        X = XA + (IX-1)*DX
        RAD = 1.D0 - (X-X0)**2/CX**2
        IF (RAD < 0.D0) CYCLE
        RY = CY * SQRT(RAD)
        RZ = CZ * SQRT(RAD)
 
        CALL EIRENE_CTELL (Y0,Z0,RY,RZ,YLIMS1,ZLIMS1,YLIMS2,ZLIMS2,RLB,
     .              PHIAN,PHIEN,IPART)
 
        DO IPRT=1,IPART
          DPH = (PHIEN(IPRT)-PHIAN(IPRT)) / FLOAT(N-1)
          Y = Y0 + RY*COS(PHIAN(IPRT))
          Z = Z0 + RZ*SIN(PHIAN(IPRT))
          CALL EIRENE_PL3D (X,Y,Z,XS,YS)
          CALL GRJMP (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          DO I = 2,N
            PHI = PHIAN(IPRT) + (I-1)*DPH
            Y = Y0 + RY*COS(PHI)
            Z = Z0 + RZ*SIN(PHI)
            CALL EIRENE_PL3D (X,Y,Z,XS,YS)
            CALL GRDRW (REAL(XS,KIND(1.E0)),REAL(YS,KIND(1.E0)))
          END DO
        END DO
      END DO ! IY
 
      CALL GRNWPN(1)
 
      RETURN
      END
