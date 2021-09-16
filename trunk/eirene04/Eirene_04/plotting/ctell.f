

      SUBROUTINE CTELL (X0,Y0,RX,RY,XLIMS1,YLIMS1,XLIMS2,YLIMS2,RLB,
     .                  PHIAN,PHIEN,IPART)

      USE PRECISION
      USE CCONA

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X0, Y0, RX, RY, XLIMS1, YLIMS1, XLIMS2,
     .                      YLIMS2, RLB
      REAL(DP), INTENT(OUT) :: PHIAN(*), PHIEN(*)
      INTEGER, INTENT(OUT) :: IPART
      REAL(DP) :: PHIANG(10)
      REAL(DP) :: PHI1, PHI2, X1, X2, Y1, Y2, XPHI, YPHI, QY1, QX2, 
     .            PHIH, PHI, QX1
      INTEGER :: IPHI, ISORT, I

      IPHI = 0

!  DETERMINE INTERSECTION POINTS OF ELLIPSE WITH BOX

!  INTERSECTION WITH X=XLIMS1
      QX1 = (XLIMS1 - X0) / (RX+EPS30)
      IF (ABS(QX1) <= 1.D0) THEN
        PHI = ACOS(QX1)
        Y1 = Y0 + RY*SIN(PHI)
        IF ((Y1 >= YLIMS1) .AND. (Y1 <= YLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI
        END IF
        PHI = PI2A - PHI
        Y2 = Y0 + RY*SIN(PHI)
        IF ((Y2 >= YLIMS1) .AND. (Y2 <= YLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI
        END IF
      END IF

!  INTERSECTION WITH X=XLIMS2
      QX2 = (XLIMS2 - X0) / (RX+EPS30)
      IF (ABS(QX2) <= 1.D0) THEN
        PHI = ACOS(QX2)
        Y1 = Y0 + RY*SIN(PHI)
        IF ((Y1 >= YLIMS1) .AND. (Y1 <= YLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI
        END IF
        PHI = PI2A - PHI
        Y2 = Y0 + RY*SIN(PHI)
        IF ((Y2 >= YLIMS1) .AND. (Y2 <= YLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI
        END IF
      END IF

!  INTERSECTION WITH Y=YLIMS1
      QY1 = (YLIMS1 - Y0) / (RY+EPS30)
      IF (ABS(QY1) <= 1.D0) THEN
        PHI1 = ASIN(QY1)
        PHI2 = PIA - PHI1
        IF (PHI1 < 0.D0) PHI1 = PHI1 + PI2A
        X1 = X0 + RX*COS(PHI1)
        IF ((X1 >= XLIMS1) .AND. (X1 <= XLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI1
        END IF
        X2 = X0 + RX*COS(PHI2)
        IF ((X2 >= XLIMS1) .AND. (X2 <= XLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI2
        END IF
      END IF

!  INTERSECTION WITH Y=YLIMS2
      QY1 = (YLIMS2 - Y0) / (RY+EPS30)
      IF (ABS(QY1) <= 1.D0) THEN
        PHI1 = ASIN(QY1)
        PHI2 = PIA - PHI1
        IF (PHI1 < 0.D0) PHI1 = PHI1 + PI2A
        X1 = X0 + RX*COS(PHI1)
        IF ((X1 >= XLIMS1) .AND. (X1 <= XLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI1
        END IF
        X2 = X0 + RX*COS(PHI2)
        IF ((X2 >= XLIMS1) .AND. (X2 <= XLIMS2)) THEN
          IPHI = IPHI + 1
          PHIANG(IPHI) = PHI2
        END IF
      END IF

      IF (IPHI == 0) THEN
        IPART = 1
        PHIAN(1) = 0.D0
        PHIEN(1) = PI2A
        RETURN
      END IF
C
C  SORT PHIANG
      ISORT = 1
      DO WHILE (ISORT > 0)
        ISORT=0
        DO I=1,IPHI-1
          IF (PHIANG(I+1).LT.PHIANG(I)) THEN
            PHIH=PHIANG(I)
            PHIANG(I)=PHIANG(I+1)
            PHIANG(I+1)=PHIH
            ISORT=ISORT+1
          ENDIF
        END DO
      END DO
C
      IPHI=IPHI+1
      PHIANG(IPHI)=PHIANG(1)
C
      IPART=0
      DO I=1,IPHI-1
        IF (ABS(PHIANG(I+1)-PHIANG(I)).LT.EPS10) CYCLE
        IF (PHIANG(I+1).LT.PHIANG(I)) PHIANG(I+1)=PHIANG(I+1)+PI2A
        PHI=0.5*(PHIANG(I)+PHIANG(I+1))
        XPHI=X0+RX*COS(PHI)
        YPHI=Y0+RY*SIN(PHI)
        IF (RLB.LT.1.5) THEN
          IF (XPHI.GE.XLIMS1.AND.XPHI.LE.XLIMS2.AND.
     .        YPHI.GE.YLIMS1.AND.YPHI.LE.YLIMS2) THEN
            IPART=IPART+1
            PHIAN(IPART)=PHIANG(I)
            PHIEN(IPART)=PHIANG(I+1)
          ENDIF
        ELSE
          IF (XPHI.LE.XLIMS1.OR.XPHI.GE.XLIMS2.OR.
     .        YPHI.LE.YLIMS1.OR.YPHI.GE.YLIMS2)THEN
            IPART=IPART+1
            PHIAN(IPART)=PHIANG(I)
            PHIEN(IPART)=PHIANG(I+1)
          ENDIF
        ENDIF
      END DO
C
      RETURN
      END
