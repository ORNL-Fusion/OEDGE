C
C
      SUBROUTINE CTQUA (A0LM,A1LM,A2LM,A3LM,A4LM,A5LM,A6LM,A7LM,A8LM,
     .                  A9LM,XLIMS1,XLIMS2,YLIMS1,YLIMS2,ZLIMS1,ZLIMS2,
     .                  RLB,RAD,CX,CY,CZ,PHIAN,PHIEN,IPART)

      USE PRECISION

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: A0LM, A1LM, A2LM, A3LM, A4LM, A5LM, A6LM,
     .                      A7LM, A8LM, A9LM, XLIMS1, XLIMS2, YLIMS1,
     .                      YLIMS2, ZLIMS1, ZLIMS2, RLB, RAD, CX, CY, CZ
      REAL(DP), INTENT(OUT) :: PHIAN(*), PHIEN(*)
      INTEGER, INTENT(OUT) :: IPART
      REAL(DP) :: PHIANG(10)
      REAL(DP) :: A0, A1, A2, A3, A4, XL1, XL2, EPS10, PI, PI180, YL1,
     .          YL2, PHI, PHIH, XPHI, YPHI, XM, YM
      INTEGER :: I, IPHI, ISORT
      LOGICAL :: LX0, LY0, LZ0, LCTX1, LCTX2, LCTY1, LCTY2
C
      EPS10=1.E-10
      PI=4.*ATAN(1.)
      PI180=PI/180.
C
      LX0=A1LM**2+A4LM**2+A7LM**2+A8LM**2.LT.EPS10
      LY0=A2LM**2+A5LM**2+A7LM**2+A9LM**2.LT.EPS10
      LZ0=A3LM**2+A6LM**2+A8LM**2+A9LM**2.LT.EPS10
C
      IF (.NOT.(LX0.OR.LY0.OR.LZ0)) THEN
C  SCHIEFER ZYLINDER
        IPART=1
        PHIAN(1)=0.
        PHIEN(1)=360.
        RETURN
      ENDIF
C
      IF ((LX0.AND.LY0).OR.(LX0.AND.LZ0).OR.(LY0.AND.LZ0)) THEN
        WRITE (6,*) ' ERROR IN CTQUA '
        WRITE (6,*) ' EQUATION DOES NOT DESCRIBE A CYLINDER '
        WRITE (6,*) ' IT DESCRIBES TWO PARALLEL PLANES '
        WRITE (6,*) ' NO PLOT IS DONE '
        IPART=0
        RETURN
      ENDIF
C
      IF (LX0) THEN
C  CYLINDER PARALLEL TO X
        A0=A0LM
        A1=A2LM
        A2=A3LM
        A3=A5LM
        A4=A6LM
        XL1=YLIMS1
        XL2=YLIMS2
        YL1=ZLIMS1
        YL2=ZLIMS2
        XM=CY
        YM=CZ
      ELSEIF (LY0) THEN
C  CYLINDER PARALLEL TO Y
        A0=A0LM
        A1=A1LM
        A2=A3LM
        A3=A4LM
        A4=A6LM
        XL1=XLIMS1
        XL2=XLIMS2
        YL1=ZLIMS1
        YL2=ZLIMS2
        XM=CX
        YM=CZ
      ELSE
C  CYLINDER PARALLEL TO Z
        A0=A0LM
        A1=A1LM
        A2=A2LM
        A3=A4LM
        A4=A5LM
        XL1=XLIMS1
        XL2=XLIMS2
        YL1=YLIMS1
        YL2=YLIMS2
        XM=CX
        YM=CY
      ENDIF
C
      LCTX1=XM-RAD.LE.XL1
      LCTX2=XM+RAD.GE.XL2
      LCTY1=YM-RAD.LE.YL1
      LCTY2=YM+RAD.GE.YL2
C
      IF (.NOT.(LCTX1.OR.LCTX2.OR.LCTY1.OR.LCTY2)) THEN
C  ZYLINDER LIEGT KOMPLETT IM QUADER
        IPART=1
        PHIAN(1)=0.
        PHIEN(1)=360.
        RETURN
      ENDIF
C
      IPHI=0
      IF (LCTX1) CALL CTCIRC (A0,A1,A2,A3,A4,XL1,YL1,YL2,XM,YM,
     .                        PHIANG,IPHI,0)
      IF (LCTX2) CALL CTCIRC (A0,A1,A2,A3,A4,XL2,YL1,YL2,XM,YM,
     .                        PHIANG,IPHI,0)
      IF (LCTY1) CALL CTCIRC (A0,A2,A1,A4,A3,YL1,XL1,XL2,YM,XM,
     .                        PHIANG,IPHI,1)
      IF (LCTY2) CALL CTCIRC (A0,A2,A1,A4,A3,YL2,XL1,XL2,YM,XM,
     .                        PHIANG,IPHI,1)
C
C  SORTIERE WINKEL
10    ISORT=0
      DO 15 I=1,IPHI-1
        IF (PHIANG(I+1).LT.PHIANG(I)) THEN
          PHIH=PHIANG(I)
          PHIANG(I)=PHIANG(I+1)
          PHIANG(I+1)=PHIH
          ISORT=ISORT+1
        ENDIF
15    CONTINUE
      IF (ISORT.GT.0) GOTO 10
C
      IPHI=IPHI+1
      PHIANG(IPHI)=PHIANG(1)
C
      IPART=0
      DO 20 I=1,IPHI-1
        IF (ABS(PHIANG(I+1)-PHIANG(I)).LT.EPS10) GOTO 20
        IF (PHIANG(I+1).LT.PHIANG(I)) PHIANG(I+1)=PHIANG(I+1)+360.
        PHI=0.5*(PHIANG(I)+PHIANG(I+1))*PI180
        XPHI=XM+RAD*COS(PHI)
        YPHI=YM+RAD*SIN(PHI)
        IF (RLB.LT.1.5) THEN
          IF (XPHI.GE.XL1.AND.XPHI.LE.XL2.AND.
     .        YPHI.GE.YL1.AND.YPHI.LE.YL2) THEN
            IPART=IPART+1
            PHIAN(IPART)=PHIANG(I)
            PHIEN(IPART)=PHIANG(I+1)
          ENDIF
        ELSE
          IF (XPHI.LE.XL1.OR.XPHI.GE.XL2.OR.YPHI.LE.YL1.OR.YPHI.GE.YL2)
     .       THEN
            IPART=IPART+1
            PHIAN(IPART)=PHIANG(I)
            PHIEN(IPART)=PHIANG(I+1)
          ENDIF
        ENDIF
20    CONTINUE
C
      RETURN
      END