C EIRENE07 COMPILATION
C ===== SOURCE: artri3.f


      FUNCTION ARTRI3 (X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3
      REAL(DP) :: A1, A2, A3, ARTRI3

      A1 = Y1*Z2 + Y3*Z1 + Y2*Z3 - Y3*Z2 - Y1*Z3 - Y2*Z1
      A2 = Z1*X2 + Z3*X1 + Z2*X3 - Z3*X2 - Z1*X3 - Z2*X1
      A3 = X1*Y2 + X3*Y1 + X2*Y3 - X3*Y2 - X1*Y3 - X2*Y1

      ARTRI3=0.5_DP * SQRT(A1**2+A2**2+A3**2)
      RETURN
      END
C ===== SOURCE: cal_vol.f


      FUNCTION CAL_VOL (P1,P2,P3,P4)
C  RETURNS VOLUME OF TETRAHEDRON DEFINED BY THE 4 POINTS P1,P2,P3,P4
      USE PRECISION
      USE CTETRA
      IMPLICIT NONE

      REAL(DP) :: CAL_VOL
      REAL(DP), INTENT(IN) :: P1(3), P2(3), P3(3), P4(3)
      REAL(DP) :: X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4
      REAL(DP) :: X1X2, X1X3, X1X4, Y1Y2, Y1Y3, Y1Y4, Z1Z2, Z1Z3, Z1Z4

      X1=P1(1)
      X2=P2(1)
      X3=P3(1)
      X4=P4(1)
      Y1=P1(2)
      Y2=P2(2)
      Y3=P3(2)
      Y4=P4(2)
      Z1=P1(3)
      Z2=P2(3)
      Z3=P3(3)
      Z4=P4(3)
      X1X2=X1-X2
      Y1Y2=Y1-Y2
      Z1Z2=Z1-Z2
      X1X3=X1-X3
      Y1Y3=Y1-Y3
      Z1Z3=Z1-Z3
      X1X4=X1-X4
      Y1Y4=Y1-Y4
      Z1Z4=Z1-Z4
      CAL_VOL = ITETHAND *
     .          (  X1X2 * Y1Y3 * Z1Z4
     .           + X1X4 * Y1Y2 * Z1Z3
     .           + X1X3 * Y1Y4 * Z1Z2
     .           - X1X4 * Y1Y3 * Z1Z2
     .           - X1X2 * Y1Y4 * Z1Z3
     .           - X1X3 * Y1Y2 * Z1Z4 ) / 6.D0

      RETURN
      END
C ===== SOURCE: insert_point.f


      SUBROUTINE INSERT_POINT (BAUM,X,Y,Z,DIST,IC)
C
      USE PRECISION
      USE PARMMOD
      USE CTETRA
      USE COMPRT, ONLY: IUNOUT
      USE module_avltree

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X,Y,Z,DIST
      INTEGER, INTENT(OUT) :: IC
      INTEGER I, ICOOR
      type(TAVLTree), pointer :: baum
      logical :: inserted

      IC = 0
      IC=NCOOR+1
      inserted=.false.
      call insert (baum, x, y, z, dist, ic, inserted)
      IF (.not.inserted) RETURN

      NCOOR=NCOOR+1
      IF (NCOOR > NCOORD) THEN
        WRITE (iunout,*) ' ALLOWED NUMBER OF COORDINATES EXCEEDED '
        WRITE (iunout,*) ' INCREASE NCOORD '
        CALL EXIT_OWN(1)
      END IF
      XTETRA(NCOOR) = X
      YTETRA(NCOOR) = Y
      ZTETRA(NCOOR) = Z
      IC = NCOOR
      RETURN
      END



C ===== SOURCE: inside.f

      LOGICAL FUNCTION INSIDE (ITET,P)
C  DETERMINES, IF A POINT P IS INSIDE OR OUTSIDE THE TETRAHEDRON NO.
C  ITET
      USE PRECISION
      USE COMUSR
      USE CTETRA
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ITET
      REAL(DP), INTENT(IN) :: P(3)
      REAL(DP) :: PC1(3),PC2(3),PC3(3),PC4(3)
      REAL(DP) :: V1,V2,V3,V4,CAL_VOL

      PC1(1:3)= (/ XTETRA(NTECK(1,ITET)), YTETRA(NTECK(1,ITET)),
     .             ZTETRA(NTECK(1,ITET)) /)
      PC2(1:3)= (/ XTETRA(NTECK(2,ITET)), YTETRA(NTECK(2,ITET)),
     .             ZTETRA(NTECK(2,ITET)) /)
      PC3(1:3)= (/ XTETRA(NTECK(3,ITET)), YTETRA(NTECK(3,ITET)),
     .             ZTETRA(NTECK(3,ITET)) /)
      PC4(1:3)= (/ XTETRA(NTECK(4,ITET)), YTETRA(NTECK(4,ITET)),
     .             ZTETRA(NTECK(4,ITET)) /)

      V1 = CAL_VOL (PC1,PC2,PC3,P)
      V2 = CAL_VOL (PC3,PC2,PC4,P)
      V3 = CAL_VOL (PC1,PC3,PC4,P)
      V4 = CAL_VOL (PC1,PC4,PC2,P)

      INSIDE = (ABS(V1+V2+V3+V4-VOL(ITET)) < 1.D-3*VOL(ITET)) .AND.
     .         (MIN(V1,V2,V3,V4) >= 0.D0)

      RETURN
      END
C ===== SOURCE: intri3.f


      LOGICAL FUNCTION INTRI3 (X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X,Y,Z)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3,
     .                        X, Y, Z
      REAL(DP) :: A1, A2, A3, A, ARTRI3
      LOGICAL :: L1, L2

      A1=ARTRI3(X1,Y1,Z1,X2,Y2,Z2,X,Y,Z)
      A2=ARTRI3(X2,Y2,Z2,X3,Y3,Z3,X,Y,Z)
      A3=ARTRI3(X3,Y3,Z3,X1,Y1,Z1,X,Y,Z)
      A=ARTRI3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)

      L1 = MIN(A1,A2,A3).GE.0.D0
      L2 = ABS(A-A1-A2-A3) < 1.D-3
      INTRI3 = L1 .AND. L2

      RETURN
      END
C ===== SOURCE: learct.f


      INTEGER FUNCTION LEARCT (X,Y,Z)

      USE PRECISION
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CTETRA

      IMPLICIT NONE

      INTEGER, PARAMETER :: NCL=30

      REAL(DP), INTENT(IN) :: X,Y,Z
      REAL(DP) :: PC1(3), PC2(3), PC3(3), PC4(3), P(3)
      REAL(DP) :: V1, V2, V3, V4, CAL_VOL
      REAL(DP), SAVE :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .                  DISTX, DISTY, DISTZ,
     .                  EPDX, EPDY, EPDZ
      REAL(DP) :: DELTAX, DELTAY, DELTAZ,
     .            XTRMIN, YTRMIN, ZTRMIN, XTRMAX, YTRMAX, ZTRMAX
      INTEGER :: ITET, IFIRST, I, J, K, IX, IY, IZ, IHEADX1, IHEADX2,
     .           IHEADY1, IHEADY2, IHEADZ1, IHEADZ2
      LOGICAL :: LG(NTET)

      TYPE :: CELL
        INTEGER :: TETNR
        TYPE(CELL),POINTER :: NEXT
      END TYPE CELL

      TYPE :: POIFELD
        TYPE (CELL),POINTER :: P
      END TYPE POIFELD

      TYPE (POIFELD) :: HELPCUR(8)
      TYPE (POIFELD),ALLOCATABLE,SAVE :: HEADS(:,:,:)
      TYPE (CELL),POINTER :: CUR

      DATA IFIRST /0/

      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
        ALLOCATE(HEADS(NCL,NCL,NCL))
        DO I=1,NCL
          DO J=1,NCL
            DO K=1,NCL
              NULLIFY(HEADS(I,J,K)%P)
            ENDDO
          ENDDO
        ENDDO
        XMIN=MINVAL(XTETRA(1:NCOOR))
        YMIN=MINVAL(YTETRA(1:NCOOR))
        ZMIN=MINVAL(ZTETRA(1:NCOOR))
        XMAX=MAXVAL(XTETRA(1:NCOOR))
        YMAX=MAXVAL(YTETRA(1:NCOOR))
        ZMAX=MAXVAL(ZTETRA(1:NCOOR))
        XMIN = XMIN *(1.-EPS5)
        XMAX = XMAX *(1.+EPS5)
        YMIN = YMIN *(1.-EPS5)
        YMAX = YMAX *(1.+EPS5)
        ZMIN = ZMIN *(1.-EPS5)
        ZMAX = ZMAX *(1.+EPS5)
        DISTX=(XMAX-XMIN)/REAL(NCL,KIND=DP)
        DISTY=(YMAX-YMIN)/REAL(NCL,KIND=DP)
        DISTZ=(ZMAX-ZMIN)/REAL(NCL,KIND=DP)
        EPDX=DISTX*EPS5
        EPDY=DISTY*EPS5
        EPDZ=DISTZ*EPS5
        DO I=1,NTET
          IF (VOL(I) < EPS30) CYCLE
          XTRMIN = MIN(XTETRA(NTECK(1,I)),XTETRA(NTECK(2,I)),
     .                 XTETRA(NTECK(3,I)),XTETRA(NTECK(4,I)))
          XTRMAX = MAX(XTETRA(NTECK(1,I)),XTETRA(NTECK(2,I)),
     .                 XTETRA(NTECK(3,I)),XTETRA(NTECK(4,I)))
          YTRMIN = MIN(YTETRA(NTECK(1,I)),YTETRA(NTECK(2,I)),
     .                 YTETRA(NTECK(3,I)),YTETRA(NTECK(4,I)))
          YTRMAX = MAX(YTETRA(NTECK(1,I)),YTETRA(NTECK(2,I)),
     .                 YTETRA(NTECK(3,I)),YTETRA(NTECK(4,I)))
          ZTRMIN = MIN(ZTETRA(NTECK(1,I)),ZTETRA(NTECK(2,I)),
     .                 ZTETRA(NTECK(3,I)),ZTETRA(NTECK(4,I)))
          ZTRMAX = MAX(ZTETRA(NTECK(1,I)),ZTETRA(NTECK(2,I)),
     .                 ZTETRA(NTECK(3,I)),ZTETRA(NTECK(4,I)))
          DELTAX=XTRMIN-XMIN
          IHEADX1=INT(DELTAX/DISTX)+1
          DELTAX=XTRMAX-XMIN
          IHEADX2=MAX(IHEADX1,INT(DELTAX/DISTX)+1)
          DELTAY=YTRMIN-YMIN
          IHEADY1=INT(DELTAY/DISTY)+1
          DELTAY=YTRMAX-YMIN
          IHEADY2=MAX(IHEADY1,INT(DELTAY/DISTY)+1)
          DELTAZ=ZTRMIN-ZMIN
          IHEADZ1=INT(DELTAZ/DISTZ)+1
          DELTAZ=ZTRMAX-ZMIN
          IHEADZ2=MAX(IHEADZ1,INT(DELTAZ/DISTZ)+1)
          DO IX=IHEADX1,IHEADX2
            DO IY=IHEADY1,IHEADY2
              DO IZ=IHEADZ1,IHEADZ2
                ALLOCATE(CUR)
                CUR%TETNR = I
                CUR%NEXT => HEADS(IX,IY,IZ)%P
                HEADS(IX,IY,IZ)%P => CUR
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      DELTAX=X-XMIN
      IHEADX2 = 0
      IF (ABS(MOD(DELTAX,DISTX)) .LT. EPDX) THEN
        IHEADX2=INT(DELTAX/DISTX)
      ENDIF
      IHEADX1=INT(DELTAX/DISTX)+1

      DELTAY=Y-YMIN
      IHEADY2 = 0
      IF (ABS(MOD(DELTAY,DISTY)) .LT. EPDY) THEN
        IHEADY2=INT(DELTAY/DISTY)
      ENDIF
      IHEADY1=INT(DELTAY/DISTY)+1

      DELTAZ=Z-ZMIN
      IHEADZ2 = 0
      IF (ABS(MOD(DELTAZ,DISTZ)) .LT. EPDZ) THEN
        IHEADZ2=INT(DELTAZ/DISTZ)
      ENDIF
      IHEADZ1=INT(DELTAZ/DISTZ)+1

      HELPCUR(1)%P => HEADS(IHEADX1,IHEADY1,IHEADZ1)%P
      IF (IHEADX2 .GT. 0) THEN
        HELPCUR(2)%P => HEADS(IHEADX2,IHEADY1,IHEADZ1)%P
      ELSE
        NULLIFY(HELPCUR(2)%P)
      ENDIF
      IF (IHEADY2 .GT. 0) THEN
        HELPCUR(3)%P => HEADS(IHEADX1,IHEADY2,IHEADZ1)%P
      ELSE
        NULLIFY(HELPCUR(3)%P)
      ENDIF
      IF ((IHEADX2 .GT. 0) .AND. (IHEADY2 .GT. 0)) THEN
        HELPCUR(4)%P => HEADS(IHEADX2,IHEADY2,IHEADZ1)%P
      ELSE
        NULLIFY(HELPCUR(4)%P)
      ENDIF

      IF (IHEADZ2 .GT. 0) THEN
        HELPCUR(5)%P => HEADS(IHEADX1,IHEADY1,IHEADZ2)%P
        IF (IHEADX2 .GT. 0) THEN
          HELPCUR(6)%P => HEADS(IHEADX2,IHEADY1,IHEADZ2)%P
        ELSE
          NULLIFY(HELPCUR(6)%P)
        ENDIF
        IF (IHEADY2 .GT. 0) THEN
          HELPCUR(7)%P => HEADS(IHEADX1,IHEADY2,IHEADZ2)%P
        ELSE
          NULLIFY(HELPCUR(7)%P)
        ENDIF
        IF ((IHEADX2 .GT. 0) .AND. (IHEADY2 .GT. 0)) THEN
          HELPCUR(8)%P => HEADS(IHEADX2,IHEADY2,IHEADZ2)%P
        ELSE
          NULLIFY(HELPCUR(8)%P)
        ENDIF
      END IF


      P(1:3) = (/ X, Y, Z /)
      LG = .FALSE.
      DO J=1,8
        DO WHILE (ASSOCIATED(HELPCUR(J)%P))
          ITET = HELPCUR(J)%P%TETNR
C  CELL I ALREADY TESTED BEFORE ?
!pb          IF (.NOT.LG(ITET) .AND. (VOL(ITET) > 1.E-6)) THEN
          IF (.NOT.LG(ITET)) THEN
            PC1(1:3)= (/ XTETRA(NTECK(1,ITET)), YTETRA(NTECK(1,ITET)),
     .                   ZTETRA(NTECK(1,ITET)) /)
            PC2(1:3)= (/ XTETRA(NTECK(2,ITET)), YTETRA(NTECK(2,ITET)),
     .                   ZTETRA(NTECK(2,ITET)) /)
            PC3(1:3)= (/ XTETRA(NTECK(3,ITET)), YTETRA(NTECK(3,ITET)),
     .                   ZTETRA(NTECK(3,ITET)) /)
            PC4(1:3)= (/ XTETRA(NTECK(4,ITET)), YTETRA(NTECK(4,ITET)),
     .                   ZTETRA(NTECK(4,ITET)) /)

            IF ((MIN(PC1(1),PC2(1),PC3(1),PC4(1)) <= X) .AND.
     .          (MAX(PC1(1),PC2(1),PC3(1),PC4(1)) >= X) .AND.
     .          (MIN(PC1(2),PC2(2),PC3(2),PC4(2)) <= Y) .AND.
     .          (MAX(PC1(2),PC2(2),PC3(2),PC4(2)) >= Y) .AND.
     .          (MIN(PC1(3),PC2(3),PC3(3),PC4(3)) <= Z) .AND.
     .          (MAX(PC1(3),PC2(3),PC3(3),PC4(3)) >= Z)) THEN

              V1 = CAL_VOL (PC1,PC2,PC3,P)
              V2 = CAL_VOL (PC3,PC2,PC4,P)
              V3 = CAL_VOL (PC1,PC3,PC4,P)
              V4 = CAL_VOL (PC1,PC4,PC2,P)

              IF ((ABS(V1+V2+V3+V4-VOL(ITET)) < 1.D-3*VOL(ITET)) .AND.
     .            (MIN(V1,V2,V3,V4) >= -EPS5*VOL(ITET))) THEN
                LEARCT=ITET
                RETURN
              END IF
            END IF
          END IF
          HELPCUR(J)%P => HELPCUR(J)%P%NEXT
        END DO
      END DO

      WRITE (iunout,*) ' POINT ',X,Y,Z
      WRITE (iunout,*) ' OUTSIDE OF ALL TETRAHEDRONS '
      LEARCT=0
      RETURN
      END
C ===== SOURCE: make_tetra.f


      SUBROUTINE MAKE_TETRA (IC1,IC2,IC3,IC4,IC5,IC6,IC7,IC8,IR,IP,IT)

      USE PRECISION
      USE PARMMOD
      USE CTETRA
      USE COMPRT, ONLY: IUNOUT
      USE CLGIN

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: IC1,IC2,IC3,IC4,IC5,IC6,IC7,IC8,IR,IP,IT

      IF (NTET+6 > NTETRA) THEN
        WRITE (iunout,*) ' ALLOWED NUMBER OF TETRAHEDRONS EXCEEDED '
        WRITE (iunout,*) ' INCREASE NTETRA '
        CALL EXIT_OWN(1)
      END IF

!      NTECK(1:4,NTET+1) = (/ IC1, IC4, IC3, IC7 /)
!      NTECK(1:4,NTET+2) = (/ IC1, IC8, IC4, IC7 /)
!      NTECK(1:4,NTET+3) = (/ IC1, IC5, IC8, IC7 /)
!      NTECK(1:4,NTET+4) = (/ IC1, IC6, IC5, IC7 /)
!      NTECK(1:4,NTET+5) = (/ IC1, IC2, IC6, IC7 /)
!      NTECK(1:4,NTET+6) = (/ IC1, IC3, IC2, IC7 /)

      NTECK(1,NTET+1) = IC1
      NTECK(2,NTET+1) = IC4
      NTECK(3,NTET+1) = IC3
      NTECK(4,NTET+1) = IC7

      NTECK(1,NTET+2) = IC1
      NTECK(2,NTET+2) = IC8
      NTECK(3,NTET+2) = IC4
      NTECK(4,NTET+2) = IC7

      NTECK(1,NTET+3) = IC1
      NTECK(2,NTET+3) = IC5
      NTECK(3,NTET+3) = IC8
      NTECK(4,NTET+3) = IC7

      NTECK(1,NTET+4) = IC1
      NTECK(2,NTET+4) = IC6
      NTECK(3,NTET+4) = IC5
      NTECK(4,NTET+4) = IC7

      NTECK(1,NTET+5) = IC1
      NTECK(2,NTET+5) = IC2
      NTECK(3,NTET+5) = IC6
      NTECK(4,NTET+5) = IC7

      NTECK(1,NTET+6) = IC1
      NTECK(2,NTET+6) = IC3
      NTECK(3,NTET+6) = IC2
      NTECK(4,NTET+6) = IC7

      CALL EINFUEGEN (IC1,NTET+1)
      CALL EINFUEGEN (IC4,NTET+1)
      CALL EINFUEGEN (IC3,NTET+1)
      CALL EINFUEGEN (IC7,NTET+1)

      CALL EINFUEGEN (IC1,NTET+2)
      CALL EINFUEGEN (IC8,NTET+2)
      CALL EINFUEGEN (IC4,NTET+2)
      CALL EINFUEGEN (IC7,NTET+2)

      CALL EINFUEGEN (IC1,NTET+3)
      CALL EINFUEGEN (IC5,NTET+3)
      CALL EINFUEGEN (IC8,NTET+3)
      CALL EINFUEGEN (IC7,NTET+3)

      CALL EINFUEGEN (IC1,NTET+4)
      CALL EINFUEGEN (IC6,NTET+4)
      CALL EINFUEGEN (IC5,NTET+4)
      CALL EINFUEGEN (IC7,NTET+4)

      CALL EINFUEGEN (IC1,NTET+5)
      CALL EINFUEGEN (IC2,NTET+5)
      CALL EINFUEGEN (IC6,NTET+5)
      CALL EINFUEGEN (IC7,NTET+5)

      CALL EINFUEGEN (IC1,NTET+6)
      CALL EINFUEGEN (IC3,NTET+6)
      CALL EINFUEGEN (IC2,NTET+6)
      CALL EINFUEGEN (IC7,NTET+6)

      IF ((IC1 == IC4) .AND. (IC5 == IC8)) THEN
        NTBAR(1:4,NTET+1:NTET+3) = -1
      ELSE
        NTBAR(2,NTET+1) = NTET+2
        NTSEITE(2,NTET+1) = 4
        NTBAR(4,NTET+1) = NTET+6
        NTSEITE(4,NTET+1) = 2

        NTBAR(2,NTET+2) = NTET+3
        NTSEITE(2,NTET+2) = 4
        NTBAR(4,NTET+2) = NTET+1
        NTSEITE(4,NTET+2) = 2

        NTBAR(2,NTET+3) = NTET+4
        NTSEITE(2,NTET+3) = 4
        NTBAR(4,NTET+3) = NTET+2
        NTSEITE(4,NTET+3) = 2

        NTBAR(4,NTET+4) = NTET+3
        NTSEITE(4,NTET+4) = 2

        NTBAR(2,NTET+6) = NTET+1
        NTSEITE(2,NTET+6) = 4
      END IF

      NTBAR(2,NTET+4) = NTET+5
      NTSEITE(2,NTET+4) = 4

      NTBAR(2,NTET+5) = NTET+6
      NTSEITE(2,NTET+5) = 4
      NTBAR(4,NTET+5) = NTET+4
      NTSEITE(4,NTET+5) = 2

      NTBAR(4,NTET+6) = NTET+5
      NTSEITE(4,NTET+6) = 2

      INMTIT(1:4,NTET+1:NTET+6) = 0

      INMTIT(1,NTET+1) = INMP3I(IR,IP,IT)
      INMTIT(3,NTET+1) = INMP2I(IR,IP+1,IT)

      INMTIT(1,NTET+2) = INMP1I(IR,IP,IT)
      INMTIT(3,NTET+2) = INMP2I(IR,IP+1,IT)

      INMTIT(1,NTET+3) = INMP1I(IR,IP,IT)
      INMTIT(3,NTET+3) = INMP3I(IR,IP,IT+1)

      INMTIT(1,NTET+4) = INMP2I(IR,IP,IT)
      INMTIT(3,NTET+4) = INMP3I(IR,IP,IT+1)

      INMTIT(1,NTET+5) = INMP2I(IR,IP,IT)
      INMTIT(3,NTET+5) = INMP1I(IR+1,IP,IT)

      INMTIT(1,NTET+6) = INMP3I(IR,IP,IT)
      INMTIT(3,NTET+6) = INMP1I(IR+1,IP,IT)

      NTET = NTET+6

      RETURN

      CONTAINS

      SUBROUTINE EINFUEGEN (IC,ITET)
        INTEGER, INTENT(IN) :: IC, ITET
        TYPE(TET_ELEM), POINTER :: CUR

        ALLOCATE (CUR)
        CUR%NOTET = ITET
        CUR%NEXT_TET => COORTET(IC)%PTET
        COORTET(IC)%PTET => CUR
        MCLSTR = MCLSTR+1
      END SUBROUTINE EINFUEGEN

      END
C ===== SOURCE: make_tetra_48.f


      SUBROUTINE MAKE_TETRA_48 (INDCO)

      USE PRECISION
      USE PARMMOD
      USE CTETRA
      USE COMPRT, ONLY: IUNOUT
      USE CLGIN

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: INDCO(27)
      INTEGER :: ITET, JTET, IS, JS
      REAL(DP), SAVE :: FEM_KOOR(3,27) = reshape(
     . (/-1._dp, -1._dp, -1._dp,
     .    0._dp, -1._dp, -1._dp,
     .    1._dp, -1._dp, -1._dp,
     .   -1._dp,  0._dp, -1._dp,
     .    0._dp,  0._dp, -1._dp,
     .    1._dp,  0._dp, -1._dp,
     .   -1._dp,  1._dp, -1._dp,
     .    0._dp,  1._dp, -1._dp,
     .    1._dp,  1._dp, -1._dp,
     .   -1._dp, -1._dp,  0._dp,
     .    0._dp, -1._dp,  0._dp,
     .    1._dp, -1._dp,  0._dp,
     .   -1._dp,  0._dp,  0._dp,
     .    0._dp,  0._dp,  0._dp,
     .    1._dp,  0._dp,  0._dp,
     .   -1._dp,  1._dp,  0._dp,
     .    0._dp,  1._dp,  0._dp,
     .    1._dp,  1._dp,  0._dp,
     .   -1._dp, -1._dp,  1._dp,
     .    0._dp, -1._dp,  1._dp,
     .    1._dp, -1._dp,  1._dp,
     .   -1._dp,  0._dp,  1._dp,
     .    0._dp,  0._dp,  1._dp,
     .    1._dp,  0._dp,  1._dp,
     .   -1._dp,  1._dp,  1._dp,
     .    0._dp,  1._dp,  1._dp,
     .    1._dp,  1._dp,  1._dp /), (/ 3, 27 /))

      IF (NTET+48 > NTETRA) THEN
        WRITE (iunout,*) ' ALLOWED NUMBER OF TETRAHEDRONS EXCEEDED '
        WRITE (iunout,*) ' INCREASE NTETRA '
        CALL EXIT_OWN(1)
      END IF

      NTECK(1,NTET+1) = INDCO(1)
      NTECK(2,NTET+1) = INDCO(2)
      NTECK(3,NTET+1) = INDCO(11)
      NTECK(4,NTET+1) = INDCO(14)
      CALL EINFUEGEN (INDCO(1),NTET+1)
      CALL EINFUEGEN (INDCO(2),NTET+1)
      CALL EINFUEGEN (INDCO(11),NTET+1)
      CALL EINFUEGEN (INDCO(14),NTET+1)
      RTCEN(NTET+1) = 0.25_DP* (FEM_KOOR(1,1)+FEM_KOOR(1,2)+
     .                          FEM_KOOR(1,11)+FEM_KOOR(1,14))
      STCEN(NTET+1) = 0.25_DP* (FEM_KOOR(2,1)+FEM_KOOR(2,2)+
     .                          FEM_KOOR(2,11)+FEM_KOOR(2,14))
      TTCEN(NTET+1) = 0.25_DP* (FEM_KOOR(3,1)+FEM_KOOR(3,2)+
     .                          FEM_KOOR(3,11)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(1),INDCO(2),INDCO(11),INDCO(14))) THEN
        NTBAR(2,NTET+1) = NTET+34
        NTSEITE(2,NTET+1) = 2
        NTBAR(3,NTET+1) = NTET+2
        NTSEITE(3,NTET+1) = 4
        NTBAR(4,NTET+1) = NTET+8
        NTSEITE(4,NTET+1) = 3
      ELSE
        NTBAR(1:4,NTET+1) = -1
        NTSEITE(1:4,NTET+1) = -1
      END IF


      NTECK(1,NTET+2) = INDCO(2)
      NTECK(2,NTET+2) = INDCO(3)
      NTECK(3,NTET+2) = INDCO(11)
      NTECK(4,NTET+2) = INDCO(14)
      CALL EINFUEGEN (INDCO(2),NTET+2)
      CALL EINFUEGEN (INDCO(3),NTET+2)
      CALL EINFUEGEN (INDCO(11),NTET+2)
      CALL EINFUEGEN (INDCO(14),NTET+2)
      RTCEN(NTET+2) = 0.25_DP* (FEM_KOOR(1,2)+FEM_KOOR(1,3)+
     .                          FEM_KOOR(1,11)+FEM_KOOR(1,14))
      STCEN(NTET+2) = 0.25_DP* (FEM_KOOR(2,2)+FEM_KOOR(2,3)+
     .                          FEM_KOOR(2,11)+FEM_KOOR(2,14))
      TTCEN(NTET+2) = 0.25_DP* (FEM_KOOR(3,2)+FEM_KOOR(3,3)+
     .                          FEM_KOOR(3,11)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(2),INDCO(3),INDCO(11),INDCO(14))) THEN
        NTBAR(2,NTET+2) = NTET+33
        NTSEITE(2,NTET+2) = 2
        NTBAR(3,NTET+2) = NTET+3
        NTSEITE(3,NTET+2) = 4
        NTBAR(4,NTET+2) = NTET+1
        NTSEITE(4,NTET+2) = 3
      ELSE
        NTBAR(1:4,NTET+2) = -1
        NTSEITE(1:4,NTET+2) = -1
      END IF


      NTECK(1,NTET+3) = INDCO(3)
      NTECK(2,NTET+3) = INDCO(12)
      NTECK(3,NTET+3) = INDCO(11)
      NTECK(4,NTET+3) = INDCO(14)
      CALL EINFUEGEN (INDCO(3),NTET+3)
      CALL EINFUEGEN (INDCO(12),NTET+3)
      CALL EINFUEGEN (INDCO(11),NTET+3)
      CALL EINFUEGEN (INDCO(14),NTET+3)
      RTCEN(NTET+3) = 0.25_DP* (FEM_KOOR(1,3)+FEM_KOOR(1,12)+
     .                          FEM_KOOR(1,11)+FEM_KOOR(1,14))
      STCEN(NTET+3) = 0.25_DP* (FEM_KOOR(2,3)+FEM_KOOR(2,12)+
     .                          FEM_KOOR(2,11)+FEM_KOOR(2,14))
      TTCEN(NTET+3) = 0.25_DP* (FEM_KOOR(3,3)+FEM_KOOR(3,12)+
     .                          FEM_KOOR(3,11)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(3),INDCO(12),INDCO(11),INDCO(14))) THEN
        NTBAR(2,NTET+3) = NTET+16
        NTSEITE(2,NTET+3) = 2
        NTBAR(3,NTET+3) = NTET+4
        NTSEITE(3,NTET+3) = 4
        NTBAR(4,NTET+3) = NTET+2
        NTSEITE(4,NTET+3) = 3
      ELSE
        NTBAR(1:4,NTET+3) = -1
        NTSEITE(1:4,NTET+3) = -1
      END IF


      NTECK(1,NTET+4) = INDCO(12)
      NTECK(2,NTET+4) = INDCO(21)
      NTECK(3,NTET+4) = INDCO(11)
      NTECK(4,NTET+4) = INDCO(14)
      CALL EINFUEGEN (INDCO(12),NTET+4)
      CALL EINFUEGEN (INDCO(21),NTET+4)
      CALL EINFUEGEN (INDCO(11),NTET+4)
      CALL EINFUEGEN (INDCO(14),NTET+4)
      RTCEN(NTET+4) = 0.25_DP* (FEM_KOOR(1,12)+FEM_KOOR(1,21)+
     .                          FEM_KOOR(1,11)+FEM_KOOR(1,14))
      STCEN(NTET+4) = 0.25_DP* (FEM_KOOR(2,12)+FEM_KOOR(2,21)+
     .                          FEM_KOOR(2,11)+FEM_KOOR(2,14))
      TTCEN(NTET+4) = 0.25_DP* (FEM_KOOR(3,12)+FEM_KOOR(3,21)+
     .                          FEM_KOOR(3,11)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(12),INDCO(21),INDCO(11),INDCO(14))) THEN
        NTBAR(2,NTET+4) = NTET+15
        NTSEITE(2,NTET+4) = 2
        NTBAR(3,NTET+4) = NTET+5
        NTSEITE(3,NTET+4) = 4
        NTBAR(4,NTET+4) = NTET+3
        NTSEITE(4,NTET+4) = 3
      ELSE
        NTBAR(1:4,NTET+4) = -1
        NTSEITE(1:4,NTET+4) = -1
      END IF


      NTECK(1,NTET+5) = INDCO(21)
      NTECK(2,NTET+5) = INDCO(20)
      NTECK(3,NTET+5) = INDCO(11)
      NTECK(4,NTET+5) = INDCO(14)
      CALL EINFUEGEN (INDCO(21),NTET+5)
      CALL EINFUEGEN (INDCO(20),NTET+5)
      CALL EINFUEGEN (INDCO(11),NTET+5)
      CALL EINFUEGEN (INDCO(14),NTET+5)
      RTCEN(NTET+5) = 0.25_DP* (FEM_KOOR(1,21)+FEM_KOOR(1,20)+
     .                          FEM_KOOR(1,11)+FEM_KOOR(1,14))
      STCEN(NTET+5) = 0.25_DP* (FEM_KOOR(2,21)+FEM_KOOR(2,20)+
     .                          FEM_KOOR(2,11)+FEM_KOOR(2,14))
      TTCEN(NTET+5) = 0.25_DP* (FEM_KOOR(3,21)+FEM_KOOR(3,20)+
     .                          FEM_KOOR(3,11)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(21),INDCO(20),INDCO(11),INDCO(14))) THEN
        NTBAR(2,NTET+5) = NTET+42
        NTSEITE(2,NTET+5) = 2
        NTBAR(3,NTET+5) = NTET+6
        NTSEITE(3,NTET+5) = 4
        NTBAR(4,NTET+5) = NTET+4
        NTSEITE(4,NTET+5) = 3
      ELSE
        NTBAR(1:4,NTET+5) = -1
        NTSEITE(1:4,NTET+5) = -1
      END IF


      NTECK(1,NTET+6) = INDCO(20)
      NTECK(2,NTET+6) = INDCO(19)
      NTECK(3,NTET+6) = INDCO(11)
      NTECK(4,NTET+6) = INDCO(14)
      CALL EINFUEGEN (INDCO(20),NTET+6)
      CALL EINFUEGEN (INDCO(19),NTET+6)
      CALL EINFUEGEN (INDCO(11),NTET+6)
      CALL EINFUEGEN (INDCO(14),NTET+6)
      RTCEN(NTET+6) = 0.25_DP* (FEM_KOOR(1,20)+FEM_KOOR(1,19)+
     .                          FEM_KOOR(1,11)+FEM_KOOR(1,14))
      STCEN(NTET+6) = 0.25_DP* (FEM_KOOR(2,20)+FEM_KOOR(2,19)+
     .                          FEM_KOOR(2,11)+FEM_KOOR(2,14))
      TTCEN(NTET+6) = 0.25_DP* (FEM_KOOR(3,20)+FEM_KOOR(3,19)+
     .                          FEM_KOOR(3,11)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(20),INDCO(19),INDCO(11),INDCO(14))) THEN
        NTBAR(2,NTET+6) = NTET+41
        NTSEITE(2,NTET+6) = 2
        NTBAR(3,NTET+6) = NTET+7
        NTSEITE(3,NTET+6) = 4
        NTBAR(4,NTET+6) = NTET+5
        NTSEITE(4,NTET+6) = 3
      ELSE
        NTBAR(1:4,NTET+6) = -1
        NTSEITE(1:4,NTET+6) = -1
      END IF


      NTECK(1,NTET+7) = INDCO(19)
      NTECK(2,NTET+7) = INDCO(10)
      NTECK(3,NTET+7) = INDCO(11)
      NTECK(4,NTET+7) = INDCO(14)
      CALL EINFUEGEN (INDCO(19),NTET+7)
      CALL EINFUEGEN (INDCO(10),NTET+7)
      CALL EINFUEGEN (INDCO(11),NTET+7)
      CALL EINFUEGEN (INDCO(14),NTET+7)
      RTCEN(NTET+7) = 0.25_DP* (FEM_KOOR(1,19)+FEM_KOOR(1,10)+
     .                          FEM_KOOR(1,11)+FEM_KOOR(1,14))
      STCEN(NTET+7) = 0.25_DP* (FEM_KOOR(2,19)+FEM_KOOR(2,10)+
     .                          FEM_KOOR(2,11)+FEM_KOOR(2,14))
      TTCEN(NTET+7) = 0.25_DP* (FEM_KOOR(3,19)+FEM_KOOR(3,10)+
     .                          FEM_KOOR(3,11)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(19),INDCO(10),INDCO(11),INDCO(14))) THEN
        NTBAR(2,NTET+7) = NTET+28
        NTSEITE(2,NTET+7) = 2
        NTBAR(3,NTET+7) = NTET+8
        NTSEITE(3,NTET+7) = 4
        NTBAR(4,NTET+7) = NTET+6
        NTSEITE(4,NTET+7) = 3
      ELSE
        NTBAR(1:4,NTET+7) = -1
        NTSEITE(1:4,NTET+7) = -1
      END IF


      NTECK(1,NTET+8) = INDCO(10)
      NTECK(2,NTET+8) = INDCO(1)
      NTECK(3,NTET+8) = INDCO(11)
      NTECK(4,NTET+8) = INDCO(14)
      CALL EINFUEGEN (INDCO(10),NTET+8)
      CALL EINFUEGEN (INDCO(1),NTET+8)
      CALL EINFUEGEN (INDCO(11),NTET+8)
      CALL EINFUEGEN (INDCO(14),NTET+8)
      RTCEN(NTET+8) = 0.25_DP* (FEM_KOOR(1,10)+FEM_KOOR(1,1)+
     .                          FEM_KOOR(1,11)+FEM_KOOR(1,14))
      STCEN(NTET+8) = 0.25_DP* (FEM_KOOR(2,10)+FEM_KOOR(2,1)+
     .                          FEM_KOOR(2,11)+FEM_KOOR(2,14))
      TTCEN(NTET+8) = 0.25_DP* (FEM_KOOR(3,10)+FEM_KOOR(3,1)+
     .                          FEM_KOOR(3,11)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(10),INDCO(1),INDCO(11),INDCO(14))) THEN
        NTBAR(2,NTET+8) = NTET+27
        NTSEITE(2,NTET+8) = 2
        NTBAR(3,NTET+8) = NTET+1
        NTSEITE(3,NTET+8) = 4
        NTBAR(4,NTET+8) = NTET+7
        NTSEITE(4,NTET+8) = 3
      ELSE
        NTBAR(1:4,NTET+8) = -1
        NTSEITE(1:4,NTET+8) = -1
      END IF


      NTECK(1,NTET+9) = INDCO(3)
      NTECK(2,NTET+9) = INDCO(6)
      NTECK(3,NTET+9) = INDCO(15)
      NTECK(4,NTET+9) = INDCO(14)
      CALL EINFUEGEN (INDCO(3),NTET+9)
      CALL EINFUEGEN (INDCO(6),NTET+9)
      CALL EINFUEGEN (INDCO(15),NTET+9)
      CALL EINFUEGEN (INDCO(14),NTET+9)
      RTCEN(NTET+9) = 0.25_DP* (FEM_KOOR(1,3)+FEM_KOOR(1,6)+
     .                          FEM_KOOR(1,15)+FEM_KOOR(1,14))
      STCEN(NTET+9) = 0.25_DP* (FEM_KOOR(2,3)+FEM_KOOR(2,6)+
     .                          FEM_KOOR(2,15)+FEM_KOOR(2,14))
      TTCEN(NTET+9) = 0.25_DP* (FEM_KOOR(3,3)+FEM_KOOR(3,6)+
     .                          FEM_KOOR(3,15)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(3),INDCO(6),INDCO(15),INDCO(14))) THEN
        NTBAR(2,NTET+9) = NTET+40
        NTSEITE(2,NTET+9) = 2
        NTBAR(3,NTET+9) = NTET+10
        NTSEITE(3,NTET+9) = 4
        NTBAR(4,NTET+9) = NTET+16
        NTSEITE(4,NTET+9) = 3
      ELSE
        NTBAR(1:4,NTET+9) = -1
        NTSEITE(1:4,NTET+9) = -1
      END IF


      NTECK(1,NTET+10) = INDCO(6)
      NTECK(2,NTET+10) = INDCO(9)
      NTECK(3,NTET+10) = INDCO(15)
      NTECK(4,NTET+10) = INDCO(14)
      CALL EINFUEGEN (INDCO(6),NTET+10)
      CALL EINFUEGEN (INDCO(9),NTET+10)
      CALL EINFUEGEN (INDCO(15),NTET+10)
      CALL EINFUEGEN (INDCO(14),NTET+10)
      RTCEN(NTET+10) = 0.25_DP* (FEM_KOOR(1,6)+FEM_KOOR(1,9)+
     .                          FEM_KOOR(1,15)+FEM_KOOR(1,14))
      STCEN(NTET+10) = 0.25_DP* (FEM_KOOR(2,6)+FEM_KOOR(2,9)+
     .                          FEM_KOOR(2,15)+FEM_KOOR(2,14))
      TTCEN(NTET+10) = 0.25_DP* (FEM_KOOR(3,6)+FEM_KOOR(3,9)+
     .                          FEM_KOOR(3,15)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(6),INDCO(9),INDCO(15),INDCO(14))) THEN
        NTBAR(2,NTET+10) = NTET+39
        NTSEITE(2,NTET+10) = 2
        NTBAR(3,NTET+10) = NTET+11
        NTSEITE(3,NTET+10) = 4
        NTBAR(4,NTET+10) = NTET+9
        NTSEITE(4,NTET+10) = 3
      ELSE
        NTBAR(1:4,NTET+10) = -1
        NTSEITE(1:4,NTET+10) = -1
      END IF


      NTECK(1,NTET+11) = INDCO(9)
      NTECK(2,NTET+11) = INDCO(18)
      NTECK(3,NTET+11) = INDCO(15)
      NTECK(4,NTET+11) = INDCO(14)
      CALL EINFUEGEN (INDCO(9),NTET+11)
      CALL EINFUEGEN (INDCO(18),NTET+11)
      CALL EINFUEGEN (INDCO(15),NTET+11)
      CALL EINFUEGEN (INDCO(14),NTET+11)
      RTCEN(NTET+11) = 0.25_DP* (FEM_KOOR(1,9)+FEM_KOOR(1,18)+
     .                          FEM_KOOR(1,15)+FEM_KOOR(1,14))
      STCEN(NTET+11) = 0.25_DP* (FEM_KOOR(2,9)+FEM_KOOR(2,18)+
     .                          FEM_KOOR(2,15)+FEM_KOOR(2,14))
      TTCEN(NTET+11) = 0.25_DP* (FEM_KOOR(3,9)+FEM_KOOR(3,18)+
     .                          FEM_KOOR(3,15)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(9),INDCO(18),INDCO(15),INDCO(14))) THEN
        NTBAR(2,NTET+11) = NTET+24
        NTSEITE(2,NTET+11) = 2
        NTBAR(3,NTET+11) = NTET+12
        NTSEITE(3,NTET+11) = 4
        NTBAR(4,NTET+11) = NTET+10
        NTSEITE(4,NTET+11) = 3
      ELSE
        NTBAR(1:4,NTET+11) = -1
        NTSEITE(1:4,NTET+11) = -1
      END IF


      NTECK(1,NTET+12) = INDCO(18)
      NTECK(2,NTET+12) = INDCO(27)
      NTECK(3,NTET+12) = INDCO(15)
      NTECK(4,NTET+12) = INDCO(14)
      CALL EINFUEGEN (INDCO(18),NTET+12)
      CALL EINFUEGEN (INDCO(27),NTET+12)
      CALL EINFUEGEN (INDCO(15),NTET+12)
      CALL EINFUEGEN (INDCO(14),NTET+12)
      RTCEN(NTET+12) = 0.25_DP* (FEM_KOOR(1,18)+FEM_KOOR(1,27)+
     .                          FEM_KOOR(1,15)+FEM_KOOR(1,14))
      STCEN(NTET+12) = 0.25_DP* (FEM_KOOR(2,18)+FEM_KOOR(2,27)+
     .                          FEM_KOOR(2,15)+FEM_KOOR(2,14))
      TTCEN(NTET+12) = 0.25_DP* (FEM_KOOR(3,18)+FEM_KOOR(3,27)+
     .                          FEM_KOOR(3,15)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(18),INDCO(27),INDCO(15),INDCO(14))) THEN
        NTBAR(2,NTET+12) = NTET+23
        NTSEITE(2,NTET+12) = 2
        NTBAR(3,NTET+12) = NTET+13
        NTSEITE(3,NTET+12) = 4
        NTBAR(4,NTET+12) = NTET+11
        NTSEITE(4,NTET+12) = 3
      ELSE
        NTBAR(1:4,NTET+12) = -1
        NTSEITE(1:4,NTET+12) = -1
      END IF


      NTECK(1,NTET+13) = INDCO(27)
      NTECK(2,NTET+13) = INDCO(24)
      NTECK(3,NTET+13) = INDCO(15)
      NTECK(4,NTET+13) = INDCO(14)
      CALL EINFUEGEN (INDCO(27),NTET+13)
      CALL EINFUEGEN (INDCO(24),NTET+13)
      CALL EINFUEGEN (INDCO(15),NTET+13)
      CALL EINFUEGEN (INDCO(14),NTET+13)
      RTCEN(NTET+13) = 0.25_DP* (FEM_KOOR(1,27)+FEM_KOOR(1,24)+
     .                          FEM_KOOR(1,15)+FEM_KOOR(1,14))
      STCEN(NTET+13) = 0.25_DP* (FEM_KOOR(2,27)+FEM_KOOR(2,24)+
     .                          FEM_KOOR(2,15)+FEM_KOOR(2,14))
      TTCEN(NTET+13) = 0.25_DP* (FEM_KOOR(3,27)+FEM_KOOR(3,24)+
     .                          FEM_KOOR(3,15)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(27),INDCO(24),INDCO(15),INDCO(14))) THEN
        NTBAR(2,NTET+13) = NTET+44
        NTSEITE(2,NTET+13) = 2
        NTBAR(3,NTET+13) = NTET+14
        NTSEITE(3,NTET+13) = 4
        NTBAR(4,NTET+13) = NTET+12
        NTSEITE(4,NTET+13) = 3
      ELSE
        NTBAR(1:4,NTET+13) = -1
        NTSEITE(1:4,NTET+13) = -1
      END IF


      NTECK(1,NTET+14) = INDCO(24)
      NTECK(2,NTET+14) = INDCO(21)
      NTECK(3,NTET+14) = INDCO(15)
      NTECK(4,NTET+14) = INDCO(14)
      CALL EINFUEGEN (INDCO(24),NTET+14)
      CALL EINFUEGEN (INDCO(21),NTET+14)
      CALL EINFUEGEN (INDCO(15),NTET+14)
      CALL EINFUEGEN (INDCO(14),NTET+14)
      RTCEN(NTET+14) = 0.25_DP* (FEM_KOOR(1,24)+FEM_KOOR(1,21)+
     .                          FEM_KOOR(1,15)+FEM_KOOR(1,14))
      STCEN(NTET+14) = 0.25_DP* (FEM_KOOR(2,24)+FEM_KOOR(2,21)+
     .                          FEM_KOOR(2,15)+FEM_KOOR(2,14))
      TTCEN(NTET+14) = 0.25_DP* (FEM_KOOR(3,24)+FEM_KOOR(3,21)+
     .                          FEM_KOOR(3,15)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(24),INDCO(21),INDCO(15),INDCO(14))) THEN
        NTBAR(2,NTET+14) = NTET+43
        NTSEITE(2,NTET+14) = 2
        NTBAR(3,NTET+14) = NTET+15
        NTSEITE(3,NTET+14) = 4
        NTBAR(4,NTET+14) = NTET+13
        NTSEITE(4,NTET+14) = 3
      ELSE
        NTBAR(1:4,NTET+14) = -1
        NTSEITE(1:4,NTET+14) = -1
      END IF


      NTECK(1,NTET+15) = INDCO(21)
      NTECK(2,NTET+15) = INDCO(12)
      NTECK(3,NTET+15) = INDCO(15)
      NTECK(4,NTET+15) = INDCO(14)
      CALL EINFUEGEN (INDCO(21),NTET+15)
      CALL EINFUEGEN (INDCO(12),NTET+15)
      CALL EINFUEGEN (INDCO(15),NTET+15)
      CALL EINFUEGEN (INDCO(14),NTET+15)
      RTCEN(NTET+15) = 0.25_DP* (FEM_KOOR(1,21)+FEM_KOOR(1,12)+
     .                          FEM_KOOR(1,15)+FEM_KOOR(1,14))
      STCEN(NTET+15) = 0.25_DP* (FEM_KOOR(2,21)+FEM_KOOR(2,12)+
     .                          FEM_KOOR(2,15)+FEM_KOOR(2,14))
      TTCEN(NTET+15) = 0.25_DP* (FEM_KOOR(3,21)+FEM_KOOR(3,12)+
     .                          FEM_KOOR(3,15)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(21),INDCO(12),INDCO(15),INDCO(14))) THEN
        NTBAR(2,NTET+15) = NTET+4
        NTSEITE(2,NTET+15) = 2
        NTBAR(3,NTET+15) = NTET+16
        NTSEITE(3,NTET+15) = 4
        NTBAR(4,NTET+15) = NTET+14
        NTSEITE(4,NTET+15) = 3
      ELSE
        NTBAR(1:4,NTET+15) = -1
        NTSEITE(1:4,NTET+15) = -1
      END IF


      NTECK(1,NTET+16) = INDCO(12)
      NTECK(2,NTET+16) = INDCO(3)
      NTECK(3,NTET+16) = INDCO(15)
      NTECK(4,NTET+16) = INDCO(14)
      CALL EINFUEGEN (INDCO(12),NTET+16)
      CALL EINFUEGEN (INDCO(3),NTET+16)
      CALL EINFUEGEN (INDCO(15),NTET+16)
      CALL EINFUEGEN (INDCO(14),NTET+16)
      RTCEN(NTET+16) = 0.25_DP* (FEM_KOOR(1,12)+FEM_KOOR(1,3)+
     .                           FEM_KOOR(1,15)+FEM_KOOR(1,14))
      STCEN(NTET+16) = 0.25_DP* (FEM_KOOR(2,12)+FEM_KOOR(2,3)+
     .                           FEM_KOOR(2,15)+FEM_KOOR(2,14))
      TTCEN(NTET+16) = 0.25_DP* (FEM_KOOR(3,12)+FEM_KOOR(3,3)+
     .                           FEM_KOOR(3,15)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(12),INDCO(3),INDCO(15),INDCO(14))) THEN
        NTBAR(2,NTET+16) = NTET+3
        NTSEITE(2,NTET+16) = 2
        NTBAR(3,NTET+16) = NTET+9
        NTSEITE(3,NTET+16) = 4
        NTBAR(4,NTET+16) = NTET+15
        NTSEITE(4,NTET+16) = 3
      ELSE
        NTBAR(1:4,NTET+16) = -1
        NTSEITE(1:4,NTET+16) = -1
      END IF


      NTECK(1,NTET+17) = INDCO(9)
      NTECK(2,NTET+17) = INDCO(8)
      NTECK(3,NTET+17) = INDCO(17)
      NTECK(4,NTET+17) = INDCO(14)
      CALL EINFUEGEN (INDCO(9),NTET+17)
      CALL EINFUEGEN (INDCO(8),NTET+17)
      CALL EINFUEGEN (INDCO(17),NTET+17)
      CALL EINFUEGEN (INDCO(14),NTET+17)
      RTCEN(NTET+17) = 0.25_DP* (FEM_KOOR(1,9)+FEM_KOOR(1,8)+
     .                           FEM_KOOR(1,17)+FEM_KOOR(1,14))
      STCEN(NTET+17) = 0.25_DP* (FEM_KOOR(2,9)+FEM_KOOR(2,8)+
     .                           FEM_KOOR(2,17)+FEM_KOOR(2,14))
      TTCEN(NTET+17) = 0.25_DP* (FEM_KOOR(3,9)+FEM_KOOR(3,8)+
     .                           FEM_KOOR(3,17)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(9),INDCO(8),INDCO(17),INDCO(14))) THEN
        NTBAR(2,NTET+17) = NTET+38
        NTSEITE(2,NTET+17) = 2
        NTBAR(3,NTET+17) = NTET+18
        NTSEITE(3,NTET+17) = 4
        NTBAR(4,NTET+17) = NTET+24
        NTSEITE(4,NTET+17) = 3
      ELSE
        NTBAR(1:4,NTET+17) = -1
        NTSEITE(1:4,NTET+17) = -1
      END IF


      NTECK(1,NTET+18) = INDCO(8)
      NTECK(2,NTET+18) = INDCO(7)
      NTECK(3,NTET+18) = INDCO(17)
      NTECK(4,NTET+18) = INDCO(14)
      CALL EINFUEGEN (INDCO(8),NTET+18)
      CALL EINFUEGEN (INDCO(7),NTET+18)
      CALL EINFUEGEN (INDCO(17),NTET+18)
      CALL EINFUEGEN (INDCO(14),NTET+18)
      RTCEN(NTET+18) = 0.25_DP* (FEM_KOOR(1,8)+FEM_KOOR(1,7)+
     .                           FEM_KOOR(1,17)+FEM_KOOR(1,14))
      STCEN(NTET+18) = 0.25_DP* (FEM_KOOR(2,8)+FEM_KOOR(2,7)+
     .                           FEM_KOOR(2,17)+FEM_KOOR(2,14))
      TTCEN(NTET+18) = 0.25_DP* (FEM_KOOR(3,8)+FEM_KOOR(3,7)+
     .                           FEM_KOOR(3,17)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(8),INDCO(7),INDCO(17),INDCO(14))) THEN
        NTBAR(2,NTET+18) = NTET+37
        NTSEITE(2,NTET+18) = 2
        NTBAR(3,NTET+18) = NTET+19
        NTSEITE(3,NTET+18) = 4
        NTBAR(4,NTET+18) = NTET+17
        NTSEITE(4,NTET+18) = 3
      ELSE
        NTBAR(1:4,NTET+18) = -1
        NTSEITE(1:4,NTET+18) = -1
      END IF


      NTECK(1,NTET+19) = INDCO(7)
      NTECK(2,NTET+19) = INDCO(16)
      NTECK(3,NTET+19) = INDCO(17)
      NTECK(4,NTET+19) = INDCO(14)
      CALL EINFUEGEN (INDCO(7),NTET+19)
      CALL EINFUEGEN (INDCO(16),NTET+19)
      CALL EINFUEGEN (INDCO(17),NTET+19)
      CALL EINFUEGEN (INDCO(14),NTET+19)
      RTCEN(NTET+19) = 0.25_DP* (FEM_KOOR(1,7)+FEM_KOOR(1,16)+
     .                           FEM_KOOR(1,17)+FEM_KOOR(1,14))
      STCEN(NTET+19) = 0.25_DP* (FEM_KOOR(2,7)+FEM_KOOR(2,16)+
     .                           FEM_KOOR(2,17)+FEM_KOOR(2,14))
      TTCEN(NTET+19) = 0.25_DP* (FEM_KOOR(3,7)+FEM_KOOR(3,16)+
     .                           FEM_KOOR(3,17)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(7),INDCO(16),INDCO(17),INDCO(14))) THEN
        NTBAR(2,NTET+19) = NTET+32
        NTSEITE(2,NTET+19) = 2
        NTBAR(3,NTET+19) = NTET+20
        NTSEITE(3,NTET+19) = 4
        NTBAR(4,NTET+19) = NTET+18
        NTSEITE(4,NTET+19) = 3
      ELSE
        NTBAR(1:4,NTET+19) = -1
        NTSEITE(1:4,NTET+19) = -1
      END IF


      NTECK(1,NTET+20) = INDCO(16)
      NTECK(2,NTET+20) = INDCO(25)
      NTECK(3,NTET+20) = INDCO(17)
      NTECK(4,NTET+20) = INDCO(14)
      CALL EINFUEGEN (INDCO(16),NTET+20)
      CALL EINFUEGEN (INDCO(25),NTET+20)
      CALL EINFUEGEN (INDCO(17),NTET+20)
      CALL EINFUEGEN (INDCO(14),NTET+20)
      RTCEN(NTET+20) = 0.25_DP* (FEM_KOOR(1,16)+FEM_KOOR(1,25)+
     .                           FEM_KOOR(1,17)+FEM_KOOR(1,14))
      STCEN(NTET+20) = 0.25_DP* (FEM_KOOR(2,16)+FEM_KOOR(2,25)+
     .                           FEM_KOOR(2,17)+FEM_KOOR(2,14))
      TTCEN(NTET+20) = 0.25_DP* (FEM_KOOR(3,16)+FEM_KOOR(3,25)+
     .                           FEM_KOOR(3,17)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(16),INDCO(25),INDCO(17),INDCO(14))) THEN
        NTBAR(2,NTET+20) = NTET+31
        NTSEITE(2,NTET+20) = 2
        NTBAR(3,NTET+20) = NTET+21
        NTSEITE(3,NTET+20) = 4
        NTBAR(4,NTET+20) = NTET+19
        NTSEITE(4,NTET+20) = 3
      ELSE
        NTBAR(1:4,NTET+20) = -1
        NTSEITE(1:4,NTET+20) = -1
      END IF


      NTECK(1,NTET+21) = INDCO(25)
      NTECK(2,NTET+21) = INDCO(26)
      NTECK(3,NTET+21) = INDCO(17)
      NTECK(4,NTET+21) = INDCO(14)
      CALL EINFUEGEN (INDCO(25),NTET+21)
      CALL EINFUEGEN (INDCO(26),NTET+21)
      CALL EINFUEGEN (INDCO(17),NTET+21)
      CALL EINFUEGEN (INDCO(14),NTET+21)
      RTCEN(NTET+21) = 0.25_DP* (FEM_KOOR(1,25)+FEM_KOOR(1,26)+
     .                           FEM_KOOR(1,17)+FEM_KOOR(1,14))
      STCEN(NTET+21) = 0.25_DP* (FEM_KOOR(2,25)+FEM_KOOR(2,26)+
     .                           FEM_KOOR(2,17)+FEM_KOOR(2,14))
      TTCEN(NTET+21) = 0.25_DP* (FEM_KOOR(3,25)+FEM_KOOR(3,26)+
     .                           FEM_KOOR(3,17)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(25),INDCO(26),INDCO(17),INDCO(14))) THEN
        NTBAR(2,NTET+21) = NTET+46
        NTSEITE(2,NTET+21) = 2
        NTBAR(3,NTET+21) = NTET+22
        NTSEITE(3,NTET+21) = 4
        NTBAR(4,NTET+21) = NTET+20
        NTSEITE(4,NTET+21) = 3
      ELSE
        NTBAR(1:4,NTET+21) = -1
        NTSEITE(1:4,NTET+21) = -1
      END IF


      NTECK(1,NTET+22) = INDCO(26)
      NTECK(2,NTET+22) = INDCO(27)
      NTECK(3,NTET+22) = INDCO(17)
      NTECK(4,NTET+22) = INDCO(14)
      CALL EINFUEGEN (INDCO(26),NTET+22)
      CALL EINFUEGEN (INDCO(27),NTET+22)
      CALL EINFUEGEN (INDCO(17),NTET+22)
      CALL EINFUEGEN (INDCO(14),NTET+22)
      RTCEN(NTET+22) = 0.25_DP* (FEM_KOOR(1,26)+FEM_KOOR(1,27)+
     .                           FEM_KOOR(1,17)+FEM_KOOR(1,14))
      STCEN(NTET+22) = 0.25_DP* (FEM_KOOR(2,26)+FEM_KOOR(2,27)+
     .                           FEM_KOOR(2,17)+FEM_KOOR(2,14))
      TTCEN(NTET+22) = 0.25_DP* (FEM_KOOR(3,26)+FEM_KOOR(3,27)+
     .                           FEM_KOOR(3,17)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(26),INDCO(27),INDCO(17),INDCO(14))) THEN
        NTBAR(2,NTET+22) = NTET+45
        NTSEITE(2,NTET+22) = 2
        NTBAR(3,NTET+22) = NTET+23
        NTSEITE(3,NTET+22) = 4
        NTBAR(4,NTET+22) = NTET+21
        NTSEITE(4,NTET+22) = 3
      ELSE
        NTBAR(1:4,NTET+22) = -1
        NTSEITE(1:4,NTET+22) = -1
      END IF


      NTECK(1,NTET+23) = INDCO(27)
      NTECK(2,NTET+23) = INDCO(18)
      NTECK(3,NTET+23) = INDCO(17)
      NTECK(4,NTET+23) = INDCO(14)
      CALL EINFUEGEN (INDCO(27),NTET+23)
      CALL EINFUEGEN (INDCO(18),NTET+23)
      CALL EINFUEGEN (INDCO(17),NTET+23)
      CALL EINFUEGEN (INDCO(14),NTET+23)
      RTCEN(NTET+23) = 0.25_DP* (FEM_KOOR(1,27)+FEM_KOOR(1,18)+
     .                           FEM_KOOR(1,17)+FEM_KOOR(1,14))
      STCEN(NTET+23) = 0.25_DP* (FEM_KOOR(2,27)+FEM_KOOR(2,18)+
     .                           FEM_KOOR(2,17)+FEM_KOOR(2,14))
      TTCEN(NTET+23) = 0.25_DP* (FEM_KOOR(3,27)+FEM_KOOR(3,18)+
     .                           FEM_KOOR(3,17)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(27),INDCO(18),INDCO(17),INDCO(14))) THEN
        NTBAR(2,NTET+23) = NTET+12
        NTSEITE(2,NTET+23) = 2
        NTBAR(3,NTET+23) = NTET+24
        NTSEITE(3,NTET+23) = 4
        NTBAR(4,NTET+23) = NTET+22
        NTSEITE(4,NTET+23) = 3
      ELSE
        NTBAR(1:4,NTET+23) = -1
        NTSEITE(1:4,NTET+23) = -1
      END IF


      NTECK(1,NTET+24) = INDCO(18)
      NTECK(2,NTET+24) = INDCO(9)
      NTECK(3,NTET+24) = INDCO(17)
      NTECK(4,NTET+24) = INDCO(14)
      CALL EINFUEGEN (INDCO(18),NTET+24)
      CALL EINFUEGEN (INDCO(9),NTET+24)
      CALL EINFUEGEN (INDCO(17),NTET+24)
      CALL EINFUEGEN (INDCO(14),NTET+24)
      RTCEN(NTET+24) = 0.25_DP* (FEM_KOOR(1,18)+FEM_KOOR(1,9)+
     .                           FEM_KOOR(1,17)+FEM_KOOR(1,14))
      STCEN(NTET+24) = 0.25_DP* (FEM_KOOR(2,18)+FEM_KOOR(2,9)+
     .                           FEM_KOOR(2,17)+FEM_KOOR(2,14))
      TTCEN(NTET+24) = 0.25_DP* (FEM_KOOR(3,18)+FEM_KOOR(3,9)+
     .                           FEM_KOOR(3,17)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(18),INDCO(9),INDCO(17),INDCO(14))) THEN
        NTBAR(2,NTET+24) = NTET+11
        NTSEITE(2,NTET+24) = 2
        NTBAR(3,NTET+24) = NTET+17
        NTSEITE(3,NTET+24) = 4
        NTBAR(4,NTET+24) = NTET+23
        NTSEITE(4,NTET+24) = 3
      ELSE
        NTBAR(1:4,NTET+24) = -1
        NTSEITE(1:4,NTET+24) = -1
      END IF


      NTECK(1,NTET+25) = INDCO(7)
      NTECK(2,NTET+25) = INDCO(4)
      NTECK(3,NTET+25) = INDCO(13)
      NTECK(4,NTET+25) = INDCO(14)
      CALL EINFUEGEN (INDCO(7),NTET+25)
      CALL EINFUEGEN (INDCO(4),NTET+25)
      CALL EINFUEGEN (INDCO(13),NTET+25)
      CALL EINFUEGEN (INDCO(14),NTET+25)
      RTCEN(NTET+25) = 0.25_DP* (FEM_KOOR(1,7)+FEM_KOOR(1,4)+
     .                           FEM_KOOR(1,13)+FEM_KOOR(1,14))
      STCEN(NTET+25) = 0.25_DP* (FEM_KOOR(2,7)+FEM_KOOR(2,4)+
     .                           FEM_KOOR(2,13)+FEM_KOOR(2,14))
      TTCEN(NTET+25) = 0.25_DP* (FEM_KOOR(3,7)+FEM_KOOR(3,4)+
     .                           FEM_KOOR(3,13)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(7),INDCO(4),INDCO(13),INDCO(14))) THEN
        NTBAR(2,NTET+25) = NTET+36
        NTSEITE(2,NTET+25) = 2
        NTBAR(3,NTET+25) = NTET+26
        NTSEITE(3,NTET+25) = 4
        NTBAR(4,NTET+25) = NTET+32
        NTSEITE(4,NTET+25) = 3
      ELSE
        NTBAR(1:4,NTET+25) = -1
        NTSEITE(1:4,NTET+25) = -1
      END IF


      NTECK(1,NTET+26) = INDCO(4)
      NTECK(2,NTET+26) = INDCO(1)
      NTECK(3,NTET+26) = INDCO(13)
      NTECK(4,NTET+26) = INDCO(14)
      CALL EINFUEGEN (INDCO(4),NTET+26)
      CALL EINFUEGEN (INDCO(1),NTET+26)
      CALL EINFUEGEN (INDCO(13),NTET+26)
      CALL EINFUEGEN (INDCO(14),NTET+26)
      RTCEN(NTET+26) = 0.25_DP* (FEM_KOOR(1,4)+FEM_KOOR(1,1)+
     .                           FEM_KOOR(1,13)+FEM_KOOR(1,14))
      STCEN(NTET+26) = 0.25_DP* (FEM_KOOR(2,4)+FEM_KOOR(2,1)+
     .                           FEM_KOOR(2,13)+FEM_KOOR(2,14))
      TTCEN(NTET+26) = 0.25_DP* (FEM_KOOR(3,4)+FEM_KOOR(3,1)+
     .                           FEM_KOOR(3,13)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(4),INDCO(1),INDCO(13),INDCO(14))) THEN
        NTBAR(2,NTET+26) = NTET+35
        NTSEITE(2,NTET+26) = 2
        NTBAR(3,NTET+26) = NTET+27
        NTSEITE(3,NTET+26) = 4
        NTBAR(4,NTET+26) = NTET+25
        NTSEITE(4,NTET+26) = 3
      ELSE
        NTBAR(1:4,NTET+26) = -1
        NTSEITE(1:4,NTET+26) = -1
      END IF


      NTECK(1,NTET+27) = INDCO(1)
      NTECK(2,NTET+27) = INDCO(10)
      NTECK(3,NTET+27) = INDCO(13)
      NTECK(4,NTET+27) = INDCO(14)
      CALL EINFUEGEN (INDCO(1),NTET+27)
      CALL EINFUEGEN (INDCO(10),NTET+27)
      CALL EINFUEGEN (INDCO(13),NTET+27)
      CALL EINFUEGEN (INDCO(14),NTET+27)
      RTCEN(NTET+27) = 0.25_DP* (FEM_KOOR(1,1)+FEM_KOOR(1,10)+
     .                           FEM_KOOR(1,13)+FEM_KOOR(1,14))
      STCEN(NTET+27) = 0.25_DP* (FEM_KOOR(2,1)+FEM_KOOR(2,10)+
     .                           FEM_KOOR(2,13)+FEM_KOOR(2,14))
      TTCEN(NTET+27) = 0.25_DP* (FEM_KOOR(3,1)+FEM_KOOR(3,10)+
     .                           FEM_KOOR(3,13)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(1),INDCO(10),INDCO(13),INDCO(14))) THEN
        NTBAR(2,NTET+27) = NTET+8
        NTSEITE(2,NTET+27) = 2
        NTBAR(3,NTET+27) = NTET+28
        NTSEITE(3,NTET+27) = 4
        NTBAR(4,NTET+27) = NTET+26
        NTSEITE(4,NTET+27) = 3
      ELSE
        NTBAR(1:4,NTET+27) = -1
        NTSEITE(1:4,NTET+27) = -1
      END IF


      NTECK(1,NTET+28) = INDCO(10)
      NTECK(2,NTET+28) = INDCO(19)
      NTECK(3,NTET+28) = INDCO(13)
      NTECK(4,NTET+28) = INDCO(14)
      CALL EINFUEGEN (INDCO(10),NTET+28)
      CALL EINFUEGEN (INDCO(19),NTET+28)
      CALL EINFUEGEN (INDCO(13),NTET+28)
      CALL EINFUEGEN (INDCO(14),NTET+28)
      RTCEN(NTET+28) = 0.25_DP* (FEM_KOOR(1,10)+FEM_KOOR(1,19)+
     .                           FEM_KOOR(1,13)+FEM_KOOR(1,14))
      STCEN(NTET+28) = 0.25_DP* (FEM_KOOR(2,10)+FEM_KOOR(2,19)+
     .                           FEM_KOOR(2,13)+FEM_KOOR(2,14))
      TTCEN(NTET+28) = 0.25_DP* (FEM_KOOR(3,10)+FEM_KOOR(3,19)+
     .                           FEM_KOOR(3,13)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(10),INDCO(19),INDCO(13),INDCO(14))) THEN
        NTBAR(2,NTET+28) = NTET+7
        NTSEITE(2,NTET+28) = 2
        NTBAR(3,NTET+28) = NTET+29
        NTSEITE(3,NTET+28) = 4
        NTBAR(4,NTET+28) = NTET+27
        NTSEITE(4,NTET+28) = 3
      ELSE
        NTBAR(1:4,NTET+28) = -1
        NTSEITE(1:4,NTET+28) = -1
      END IF

      NTECK(1,NTET+29) = INDCO(19)
      NTECK(2,NTET+29) = INDCO(22)
      NTECK(3,NTET+29) = INDCO(13)
      NTECK(4,NTET+29) = INDCO(14)
      CALL EINFUEGEN (INDCO(19),NTET+29)
      CALL EINFUEGEN (INDCO(22),NTET+29)
      CALL EINFUEGEN (INDCO(13),NTET+29)
      CALL EINFUEGEN (INDCO(14),NTET+29)
      RTCEN(NTET+29) = 0.25_DP* (FEM_KOOR(1,19)+FEM_KOOR(1,22)+
     .                           FEM_KOOR(1,13)+FEM_KOOR(1,14))
      STCEN(NTET+29) = 0.25_DP* (FEM_KOOR(2,19)+FEM_KOOR(2,22)+
     .                           FEM_KOOR(2,13)+FEM_KOOR(2,14))
      TTCEN(NTET+29) = 0.25_DP* (FEM_KOOR(3,19)+FEM_KOOR(3,22)+
     .                           FEM_KOOR(3,13)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(19),INDCO(22),INDCO(13),INDCO(14))) THEN
        NTBAR(2,NTET+29) = NTET+48
        NTSEITE(2,NTET+29) = 2
        NTBAR(3,NTET+29) = NTET+30
        NTSEITE(3,NTET+29) = 4
        NTBAR(4,NTET+29) = NTET+28
        NTSEITE(4,NTET+29) = 3
      ELSE
        NTBAR(1:4,NTET+29) = -1
        NTSEITE(1:4,NTET+29) = -1
      END IF


      NTECK(1,NTET+30) = INDCO(22)
      NTECK(2,NTET+30) = INDCO(25)
      NTECK(3,NTET+30) = INDCO(13)
      NTECK(4,NTET+30) = INDCO(14)
      CALL EINFUEGEN (INDCO(22),NTET+30)
      CALL EINFUEGEN (INDCO(25),NTET+30)
      CALL EINFUEGEN (INDCO(13),NTET+30)
      CALL EINFUEGEN (INDCO(14),NTET+30)
      RTCEN(NTET+30) = 0.25_DP* (FEM_KOOR(1,22)+FEM_KOOR(1,25)+
     .                           FEM_KOOR(1,13)+FEM_KOOR(1,14))
      STCEN(NTET+30) = 0.25_DP* (FEM_KOOR(2,22)+FEM_KOOR(2,25)+
     .                           FEM_KOOR(2,13)+FEM_KOOR(2,14))
      TTCEN(NTET+30) = 0.25_DP* (FEM_KOOR(3,22)+FEM_KOOR(3,25)+
     .                           FEM_KOOR(3,13)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(22),INDCO(25),INDCO(13),INDCO(14))) THEN
        NTBAR(2,NTET+30) = NTET+47
        NTSEITE(2,NTET+30) = 2
        NTBAR(3,NTET+30) = NTET+31
        NTSEITE(3,NTET+30) = 4
        NTBAR(4,NTET+30) = NTET+29
        NTSEITE(4,NTET+30) = 3
      ELSE
        NTBAR(1:4,NTET+30) = -1
        NTSEITE(1:4,NTET+30) = -1
      END IF


      NTECK(1,NTET+31) = INDCO(25)
      NTECK(2,NTET+31) = INDCO(16)
      NTECK(3,NTET+31) = INDCO(13)
      NTECK(4,NTET+31) = INDCO(14)
      CALL EINFUEGEN (INDCO(25),NTET+31)
      CALL EINFUEGEN (INDCO(16),NTET+31)
      CALL EINFUEGEN (INDCO(13),NTET+31)
      CALL EINFUEGEN (INDCO(14),NTET+31)
      RTCEN(NTET+31) = 0.25_DP* (FEM_KOOR(1,25)+FEM_KOOR(1,16)+
     .                           FEM_KOOR(1,13)+FEM_KOOR(1,14))
      STCEN(NTET+31) = 0.25_DP* (FEM_KOOR(2,25)+FEM_KOOR(2,16)+
     .                           FEM_KOOR(2,13)+FEM_KOOR(2,14))
      TTCEN(NTET+31) = 0.25_DP* (FEM_KOOR(3,25)+FEM_KOOR(3,16)+
     .                           FEM_KOOR(3,13)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(25),INDCO(16),INDCO(13),INDCO(14))) THEN
        NTBAR(2,NTET+31) = NTET+20
        NTSEITE(2,NTET+31) = 2
        NTBAR(3,NTET+31) = NTET+32
        NTSEITE(3,NTET+31) = 4
        NTBAR(4,NTET+31) = NTET+30
        NTSEITE(4,NTET+31) = 3
      ELSE
        NTBAR(1:4,NTET+31) = -1
        NTSEITE(1:4,NTET+31) = -1
      END IF


      NTECK(1,NTET+32) = INDCO(16)
      NTECK(2,NTET+32) = INDCO(7)
      NTECK(3,NTET+32) = INDCO(13)
      NTECK(4,NTET+32) = INDCO(14)
      CALL EINFUEGEN (INDCO(16),NTET+32)
      CALL EINFUEGEN (INDCO(7),NTET+32)
      CALL EINFUEGEN (INDCO(13),NTET+32)
      CALL EINFUEGEN (INDCO(14),NTET+32)
      RTCEN(NTET+32) = 0.25_DP* (FEM_KOOR(1,16)+FEM_KOOR(1,7)+
     .                           FEM_KOOR(1,13)+FEM_KOOR(1,14))
      STCEN(NTET+32) = 0.25_DP* (FEM_KOOR(2,16)+FEM_KOOR(2,7)+
     .                           FEM_KOOR(2,13)+FEM_KOOR(2,14))
      TTCEN(NTET+32) = 0.25_DP* (FEM_KOOR(3,16)+FEM_KOOR(3,7)+
     .                           FEM_KOOR(3,13)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(16),INDCO(7),INDCO(13),INDCO(14))) THEN
        NTBAR(2,NTET+32) = NTET+19
        NTSEITE(2,NTET+32) = 2
        NTBAR(3,NTET+32) = NTET+25
        NTSEITE(3,NTET+32) = 4
        NTBAR(4,NTET+32) = NTET+31
        NTSEITE(4,NTET+32) = 3
      ELSE
        NTBAR(1:4,NTET+32) = -1
        NTSEITE(1:4,NTET+32) = -1
      END IF


      NTECK(1,NTET+33) = INDCO(3)
      NTECK(2,NTET+33) = INDCO(2)
      NTECK(3,NTET+33) = INDCO(5)
      NTECK(4,NTET+33) = INDCO(14)
      CALL EINFUEGEN (INDCO(3),NTET+33)
      CALL EINFUEGEN (INDCO(2),NTET+33)
      CALL EINFUEGEN (INDCO(5),NTET+33)
      CALL EINFUEGEN (INDCO(14),NTET+33)
      RTCEN(NTET+33) = 0.25_DP* (FEM_KOOR(1,3)+FEM_KOOR(1,2)+
     .                           FEM_KOOR(1,5)+FEM_KOOR(1,14))
      STCEN(NTET+33) = 0.25_DP* (FEM_KOOR(2,3)+FEM_KOOR(2,2)+
     .                           FEM_KOOR(2,5)+FEM_KOOR(2,14))
      TTCEN(NTET+33) = 0.25_DP* (FEM_KOOR(3,3)+FEM_KOOR(3,2)+
     .                           FEM_KOOR(3,5)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(3),INDCO(2),INDCO(5),INDCO(14))) THEN
        NTBAR(2,NTET+33) = NTET+2
        NTSEITE(2,NTET+33) = 2
        NTBAR(3,NTET+33) = NTET+34
        NTSEITE(3,NTET+33) = 4
        NTBAR(4,NTET+33) = NTET+40
        NTSEITE(4,NTET+33) = 3
      ELSE
        NTBAR(1:4,NTET+33) = -1
        NTSEITE(1:4,NTET+33) = -1
      END IF


      NTECK(1,NTET+34) = INDCO(2)
      NTECK(2,NTET+34) = INDCO(1)
      NTECK(3,NTET+34) = INDCO(5)
      NTECK(4,NTET+34) = INDCO(14)
      CALL EINFUEGEN (INDCO(2),NTET+34)
      CALL EINFUEGEN (INDCO(1),NTET+34)
      CALL EINFUEGEN (INDCO(5),NTET+34)
      CALL EINFUEGEN (INDCO(14),NTET+34)
      RTCEN(NTET+34) = 0.25_DP* (FEM_KOOR(1,2)+FEM_KOOR(1,1)+
     .                           FEM_KOOR(1,5)+FEM_KOOR(1,14))
      STCEN(NTET+34) = 0.25_DP* (FEM_KOOR(2,2)+FEM_KOOR(2,1)+
     .                           FEM_KOOR(2,5)+FEM_KOOR(2,14))
      TTCEN(NTET+34) = 0.25_DP* (FEM_KOOR(3,2)+FEM_KOOR(3,1)+
     .                           FEM_KOOR(3,5)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(2),INDCO(1),INDCO(5),INDCO(14))) THEN
        NTBAR(2,NTET+34) = NTET+1
        NTSEITE(2,NTET+34) = 2
        NTBAR(3,NTET+34) = NTET+35
        NTSEITE(3,NTET+34) = 4
        NTBAR(4,NTET+34) = NTET+33
        NTSEITE(4,NTET+34) = 3
      ELSE
        NTBAR(1:4,NTET+34) = -1
        NTSEITE(1:4,NTET+34) = -1
      END IF


      NTECK(1,NTET+35) = INDCO(1)
      NTECK(2,NTET+35) = INDCO(4)
      NTECK(3,NTET+35) = INDCO(5)
      NTECK(4,NTET+35) = INDCO(14)
      CALL EINFUEGEN (INDCO(1),NTET+35)
      CALL EINFUEGEN (INDCO(4),NTET+35)
      CALL EINFUEGEN (INDCO(5),NTET+35)
      CALL EINFUEGEN (INDCO(14),NTET+35)
      RTCEN(NTET+35) = 0.25_DP* (FEM_KOOR(1,1)+FEM_KOOR(1,4)+
     .                           FEM_KOOR(1,5)+FEM_KOOR(1,14))
      STCEN(NTET+35) = 0.25_DP* (FEM_KOOR(2,1)+FEM_KOOR(2,4)+
     .                           FEM_KOOR(2,5)+FEM_KOOR(2,14))
      TTCEN(NTET+35) = 0.25_DP* (FEM_KOOR(3,1)+FEM_KOOR(3,4)+
     .                           FEM_KOOR(3,5)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(1),INDCO(4),INDCO(5),INDCO(14))) THEN
        NTBAR(2,NTET+35) = NTET+26
        NTSEITE(2,NTET+35) = 2
        NTBAR(3,NTET+35) = NTET+36
        NTSEITE(3,NTET+35) = 4
        NTBAR(4,NTET+35) = NTET+34
        NTSEITE(4,NTET+35) = 3
      ELSE
        NTBAR(1:4,NTET+35) = -1
        NTSEITE(1:4,NTET+35) = -1
      END IF


      NTECK(1,NTET+36) = INDCO(4)
      NTECK(2,NTET+36) = INDCO(7)
      NTECK(3,NTET+36) = INDCO(5)
      NTECK(4,NTET+36) = INDCO(14)
      CALL EINFUEGEN (INDCO(4),NTET+36)
      CALL EINFUEGEN (INDCO(7),NTET+36)
      CALL EINFUEGEN (INDCO(5),NTET+36)
      CALL EINFUEGEN (INDCO(14),NTET+36)
      RTCEN(NTET+36) = 0.25_DP* (FEM_KOOR(1,4)+FEM_KOOR(1,7)+
     .                           FEM_KOOR(1,5)+FEM_KOOR(1,14))
      STCEN(NTET+36) = 0.25_DP* (FEM_KOOR(2,4)+FEM_KOOR(2,7)+
     .                           FEM_KOOR(2,5)+FEM_KOOR(2,14))
      TTCEN(NTET+36) = 0.25_DP* (FEM_KOOR(3,4)+FEM_KOOR(3,7)+
     .                           FEM_KOOR(3,5)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(4),INDCO(7),INDCO(5),INDCO(14))) THEN
        NTBAR(2,NTET+36) = NTET+25
        NTSEITE(2,NTET+36) = 2
        NTBAR(3,NTET+36) = NTET+37
        NTSEITE(3,NTET+36) = 4
        NTBAR(4,NTET+36) = NTET+35
        NTSEITE(4,NTET+36) = 3
      ELSE
        NTBAR(1:4,NTET+36) = -1
        NTSEITE(1:4,NTET+36) = -1
      END IF


      NTECK(1,NTET+37) = INDCO(7)
      NTECK(2,NTET+37) = INDCO(8)
      NTECK(3,NTET+37) = INDCO(5)
      NTECK(4,NTET+37) = INDCO(14)
      CALL EINFUEGEN (INDCO(7),NTET+37)
      CALL EINFUEGEN (INDCO(8),NTET+37)
      CALL EINFUEGEN (INDCO(5),NTET+37)
      CALL EINFUEGEN (INDCO(14),NTET+37)
      RTCEN(NTET+37) = 0.25_DP* (FEM_KOOR(1,7)+FEM_KOOR(1,8)+
     .                           FEM_KOOR(1,5)+FEM_KOOR(1,14))
      STCEN(NTET+37) = 0.25_DP* (FEM_KOOR(2,7)+FEM_KOOR(2,8)+
     .                           FEM_KOOR(2,5)+FEM_KOOR(2,14))
      TTCEN(NTET+37) = 0.25_DP* (FEM_KOOR(3,7)+FEM_KOOR(3,8)+
     .                           FEM_KOOR(3,5)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(7),INDCO(8),INDCO(5),INDCO(14))) THEN
        NTBAR(2,NTET+37) = NTET+18
        NTSEITE(2,NTET+37) = 2
        NTBAR(3,NTET+37) = NTET+38
        NTSEITE(3,NTET+37) = 4
        NTBAR(4,NTET+37) = NTET+36
        NTSEITE(4,NTET+37) = 3
      ELSE
        NTBAR(1:4,NTET+37) = -1
        NTSEITE(1:4,NTET+37) = -1
      END IF


      NTECK(1,NTET+38) = INDCO(8)
      NTECK(2,NTET+38) = INDCO(9)
      NTECK(3,NTET+38) = INDCO(5)
      NTECK(4,NTET+38) = INDCO(14)
      CALL EINFUEGEN (INDCO(8),NTET+38)
      CALL EINFUEGEN (INDCO(9),NTET+38)
      CALL EINFUEGEN (INDCO(5),NTET+38)
      CALL EINFUEGEN (INDCO(14),NTET+38)
      RTCEN(NTET+38) = 0.25_DP* (FEM_KOOR(1,8)+FEM_KOOR(1,9)+
     .                           FEM_KOOR(1,5)+FEM_KOOR(1,14))
      STCEN(NTET+38) = 0.25_DP* (FEM_KOOR(2,8)+FEM_KOOR(2,9)+
     .                           FEM_KOOR(2,5)+FEM_KOOR(2,14))
      TTCEN(NTET+38) = 0.25_DP* (FEM_KOOR(3,8)+FEM_KOOR(3,9)+
     .                           FEM_KOOR(3,5)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(8),INDCO(9),INDCO(5),INDCO(14))) THEN
        NTBAR(2,NTET+38) = NTET+17
        NTSEITE(2,NTET+38) = 2
        NTBAR(3,NTET+38) = NTET+39
        NTSEITE(3,NTET+38) = 4
        NTBAR(4,NTET+38) = NTET+37
        NTSEITE(4,NTET+38) = 3
      ELSE
        NTBAR(1:4,NTET+38) = -1
        NTSEITE(1:4,NTET+38) = -1
      END IF


      NTECK(1,NTET+39) = INDCO(9)
      NTECK(2,NTET+39) = INDCO(6)
      NTECK(3,NTET+39) = INDCO(5)
      NTECK(4,NTET+39) = INDCO(14)
      CALL EINFUEGEN (INDCO(9),NTET+39)
      CALL EINFUEGEN (INDCO(6),NTET+39)
      CALL EINFUEGEN (INDCO(5),NTET+39)
      CALL EINFUEGEN (INDCO(14),NTET+39)
      RTCEN(NTET+39) = 0.25_DP* (FEM_KOOR(1,9)+FEM_KOOR(1,6)+
     .                           FEM_KOOR(1,5)+FEM_KOOR(1,14))
      STCEN(NTET+39) = 0.25_DP* (FEM_KOOR(2,9)+FEM_KOOR(2,6)+
     .                           FEM_KOOR(2,5)+FEM_KOOR(2,14))
      TTCEN(NTET+39) = 0.25_DP* (FEM_KOOR(3,9)+FEM_KOOR(3,6)+
     .                           FEM_KOOR(3,5)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(9),INDCO(6),INDCO(5),INDCO(14))) THEN
        NTBAR(2,NTET+39) = NTET+10
        NTSEITE(2,NTET+39) = 2
        NTBAR(3,NTET+39) = NTET+40
        NTSEITE(3,NTET+39) = 4
        NTBAR(4,NTET+39) = NTET+38
        NTSEITE(4,NTET+39) = 3
      ELSE
        NTBAR(1:4,NTET+39) = -1
        NTSEITE(1:4,NTET+39) = -1
      END IF


      NTECK(1,NTET+40) = INDCO(6)
      NTECK(2,NTET+40) = INDCO(3)
      NTECK(3,NTET+40) = INDCO(5)
      NTECK(4,NTET+40) = INDCO(14)
      CALL EINFUEGEN (INDCO(6),NTET+40)
      CALL EINFUEGEN (INDCO(3),NTET+40)
      CALL EINFUEGEN (INDCO(5),NTET+40)
      CALL EINFUEGEN (INDCO(14),NTET+40)
      RTCEN(NTET+40) = 0.25_DP* (FEM_KOOR(1,6)+FEM_KOOR(1,3)+
     .                           FEM_KOOR(1,5)+FEM_KOOR(1,14))
      STCEN(NTET+40) = 0.25_DP* (FEM_KOOR(2,6)+FEM_KOOR(2,3)+
     .                           FEM_KOOR(2,5)+FEM_KOOR(2,14))
      TTCEN(NTET+40) = 0.25_DP* (FEM_KOOR(3,6)+FEM_KOOR(3,3)+
     .                           FEM_KOOR(3,5)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(6),INDCO(3),INDCO(5),INDCO(14))) THEN
        NTBAR(2,NTET+40) = NTET+9
        NTSEITE(2,NTET+40) = 2
        NTBAR(3,NTET+40) = NTET+33
        NTSEITE(3,NTET+40) = 4
        NTBAR(4,NTET+40) = NTET+39
        NTSEITE(4,NTET+40) = 3
      ELSE
        NTBAR(1:4,NTET+40) = -1
        NTSEITE(1:4,NTET+40) = -1
      END IF


      NTECK(1,NTET+41) = INDCO(19)
      NTECK(2,NTET+41) = INDCO(20)
      NTECK(3,NTET+41) = INDCO(23)
      NTECK(4,NTET+41) = INDCO(14)
      CALL EINFUEGEN (INDCO(19),NTET+41)
      CALL EINFUEGEN (INDCO(20),NTET+41)
      CALL EINFUEGEN (INDCO(23),NTET+41)
      CALL EINFUEGEN (INDCO(14),NTET+41)
      RTCEN(NTET+41) = 0.25_DP* (FEM_KOOR(1,19)+FEM_KOOR(1,20)+
     .                           FEM_KOOR(1,23)+FEM_KOOR(1,14))
      STCEN(NTET+41) = 0.25_DP* (FEM_KOOR(2,19)+FEM_KOOR(2,20)+
     .                           FEM_KOOR(2,23)+FEM_KOOR(2,14))
      TTCEN(NTET+41) = 0.25_DP* (FEM_KOOR(3,19)+FEM_KOOR(3,20)+
     .                           FEM_KOOR(3,23)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(19),INDCO(20),INDCO(23),INDCO(14))) THEN
        NTBAR(2,NTET+41) = NTET+6
        NTSEITE(2,NTET+41) = 2
        NTBAR(3,NTET+41) = NTET+42
        NTSEITE(3,NTET+41) = 4
        NTBAR(4,NTET+41) = NTET+48
        NTSEITE(4,NTET+41) = 3
      ELSE
        NTBAR(1:4,NTET+41) = -1
        NTSEITE(1:4,NTET+41) = -1
      END IF


      NTECK(1,NTET+42) = INDCO(20)
      NTECK(2,NTET+42) = INDCO(21)
      NTECK(3,NTET+42) = INDCO(23)
      NTECK(4,NTET+42) = INDCO(14)
      CALL EINFUEGEN (INDCO(20),NTET+42)
      CALL EINFUEGEN (INDCO(21),NTET+42)
      CALL EINFUEGEN (INDCO(23),NTET+42)
      CALL EINFUEGEN (INDCO(14),NTET+42)
      RTCEN(NTET+42) = 0.25_DP* (FEM_KOOR(1,20)+FEM_KOOR(1,21)+
     .                           FEM_KOOR(1,23)+FEM_KOOR(1,14))
      STCEN(NTET+42) = 0.25_DP* (FEM_KOOR(2,20)+FEM_KOOR(2,21)+
     .                           FEM_KOOR(2,23)+FEM_KOOR(2,14))
      TTCEN(NTET+42) = 0.25_DP* (FEM_KOOR(3,20)+FEM_KOOR(3,21)+
     .                           FEM_KOOR(3,23)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(20),INDCO(21),INDCO(23),INDCO(14))) THEN
        NTBAR(2,NTET+42) = NTET+5
        NTSEITE(2,NTET+42) = 2
        NTBAR(3,NTET+42) = NTET+43
        NTSEITE(3,NTET+42) = 4
        NTBAR(4,NTET+42) = NTET+41
        NTSEITE(4,NTET+42) = 3
      ELSE
        NTBAR(1:4,NTET+42) = -1
        NTSEITE(1:4,NTET+42) = -1
      END IF


      NTECK(1,NTET+43) = INDCO(21)
      NTECK(2,NTET+43) = INDCO(24)
      NTECK(3,NTET+43) = INDCO(23)
      NTECK(4,NTET+43) = INDCO(14)
      CALL EINFUEGEN (INDCO(21),NTET+43)
      CALL EINFUEGEN (INDCO(24),NTET+43)
      CALL EINFUEGEN (INDCO(23),NTET+43)
      CALL EINFUEGEN (INDCO(14),NTET+43)
      RTCEN(NTET+43) = 0.25_DP* (FEM_KOOR(1,21)+FEM_KOOR(1,24)+
     .                           FEM_KOOR(1,23)+FEM_KOOR(1,14))
      STCEN(NTET+43) = 0.25_DP* (FEM_KOOR(2,21)+FEM_KOOR(2,24)+
     .                           FEM_KOOR(2,23)+FEM_KOOR(2,14))
      TTCEN(NTET+43) = 0.25_DP* (FEM_KOOR(3,21)+FEM_KOOR(3,24)+
     .                           FEM_KOOR(3,23)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(21),INDCO(24),INDCO(23),INDCO(14))) THEN
        NTBAR(2,NTET+43) = NTET+14
        NTSEITE(2,NTET+43) = 2
        NTBAR(3,NTET+43) = NTET+44
        NTSEITE(3,NTET+43) = 4
        NTBAR(4,NTET+43) = NTET+42
        NTSEITE(4,NTET+43) = 3
      ELSE
        NTBAR(1:4,NTET+43) = -1
        NTSEITE(1:4,NTET+43) = -1
      END IF


      NTECK(1,NTET+44) = INDCO(24)
      NTECK(2,NTET+44) = INDCO(27)
      NTECK(3,NTET+44) = INDCO(23)
      NTECK(4,NTET+44) = INDCO(14)
      CALL EINFUEGEN (INDCO(24),NTET+44)
      CALL EINFUEGEN (INDCO(27),NTET+44)
      CALL EINFUEGEN (INDCO(23),NTET+44)
      CALL EINFUEGEN (INDCO(14),NTET+44)
      RTCEN(NTET+44) = 0.25_DP* (FEM_KOOR(1,24)+FEM_KOOR(1,27)+
     .                           FEM_KOOR(1,23)+FEM_KOOR(1,14))
      STCEN(NTET+44) = 0.25_DP* (FEM_KOOR(2,24)+FEM_KOOR(2,27)+
     .                           FEM_KOOR(2,23)+FEM_KOOR(2,14))
      TTCEN(NTET+44) = 0.25_DP* (FEM_KOOR(3,24)+FEM_KOOR(3,27)+
     .                           FEM_KOOR(3,23)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(24),INDCO(27),INDCO(23),INDCO(14))) THEN
        NTBAR(2,NTET+44) = NTET+13
        NTSEITE(2,NTET+44) = 2
        NTBAR(3,NTET+44) = NTET+45
        NTSEITE(3,NTET+44) = 4
        NTBAR(4,NTET+44) = NTET+43
        NTSEITE(4,NTET+44) = 3
      ELSE
        NTBAR(1:4,NTET+44) = -1
        NTSEITE(1:4,NTET+44) = -1
      END IF


      NTECK(1,NTET+45) = INDCO(27)
      NTECK(2,NTET+45) = INDCO(26)
      NTECK(3,NTET+45) = INDCO(23)
      NTECK(4,NTET+45) = INDCO(14)
      CALL EINFUEGEN (INDCO(27),NTET+45)
      CALL EINFUEGEN (INDCO(26),NTET+45)
      CALL EINFUEGEN (INDCO(23),NTET+45)
      CALL EINFUEGEN (INDCO(14),NTET+45)
      RTCEN(NTET+45) = 0.25_DP* (FEM_KOOR(1,27)+FEM_KOOR(1,26)+
     .                           FEM_KOOR(1,23)+FEM_KOOR(1,14))
      STCEN(NTET+45) = 0.25_DP* (FEM_KOOR(2,27)+FEM_KOOR(2,26)+
     .                           FEM_KOOR(2,23)+FEM_KOOR(2,14))
      TTCEN(NTET+45) = 0.25_DP* (FEM_KOOR(3,27)+FEM_KOOR(3,26)+
     .                           FEM_KOOR(3,23)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(27),INDCO(26),INDCO(23),INDCO(14))) THEN
        NTBAR(2,NTET+45) = NTET+22
        NTSEITE(2,NTET+45) = 2
        NTBAR(3,NTET+45) = NTET+46
        NTSEITE(3,NTET+45) = 4
        NTBAR(4,NTET+45) = NTET+44
        NTSEITE(4,NTET+45) = 3
      ELSE
        NTBAR(1:4,NTET+45) = -1
        NTSEITE(1:4,NTET+45) = -1
      END IF


      NTECK(1,NTET+46) = INDCO(26)
      NTECK(2,NTET+46) = INDCO(25)
      NTECK(3,NTET+46) = INDCO(23)
      NTECK(4,NTET+46) = INDCO(14)
      CALL EINFUEGEN (INDCO(26),NTET+46)
      CALL EINFUEGEN (INDCO(25),NTET+46)
      CALL EINFUEGEN (INDCO(23),NTET+46)
      CALL EINFUEGEN (INDCO(14),NTET+46)
      RTCEN(NTET+46) = 0.25_DP* (FEM_KOOR(1,26)+FEM_KOOR(1,25)+
     .                           FEM_KOOR(1,23)+FEM_KOOR(1,14))
      STCEN(NTET+46) = 0.25_DP* (FEM_KOOR(2,26)+FEM_KOOR(2,25)+
     .                           FEM_KOOR(2,23)+FEM_KOOR(2,14))
      TTCEN(NTET+46) = 0.25_DP* (FEM_KOOR(3,26)+FEM_KOOR(3,25)+
     .                           FEM_KOOR(3,23)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(26),INDCO(25),INDCO(23),INDCO(14))) THEN
        NTBAR(2,NTET+46) = NTET+21
        NTSEITE(2,NTET+46) = 2
        NTBAR(3,NTET+46) = NTET+47
        NTSEITE(3,NTET+46) = 4
        NTBAR(4,NTET+46) = NTET+45
        NTSEITE(4,NTET+46) = 3
      ELSE
        NTBAR(1:4,NTET+46) = -1
        NTSEITE(1:4,NTET+46) = -1
      END IF


      NTECK(1,NTET+47) = INDCO(25)
      NTECK(2,NTET+47) = INDCO(22)
      NTECK(3,NTET+47) = INDCO(23)
      NTECK(4,NTET+47) = INDCO(14)
      CALL EINFUEGEN (INDCO(25),NTET+47)
      CALL EINFUEGEN (INDCO(22),NTET+47)
      CALL EINFUEGEN (INDCO(23),NTET+47)
      CALL EINFUEGEN (INDCO(14),NTET+47)
      RTCEN(NTET+47) = 0.25_DP* (FEM_KOOR(1,25)+FEM_KOOR(1,22)+
     .                           FEM_KOOR(1,23)+FEM_KOOR(1,14))
      STCEN(NTET+47) = 0.25_DP* (FEM_KOOR(2,25)+FEM_KOOR(2,22)+
     .                           FEM_KOOR(2,23)+FEM_KOOR(2,14))
      TTCEN(NTET+47) = 0.25_DP* (FEM_KOOR(3,25)+FEM_KOOR(3,22)+
     .                           FEM_KOOR(3,23)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(25),INDCO(22),INDCO(23),INDCO(14))) THEN
        NTBAR(2,NTET+47) = NTET+30
        NTSEITE(2,NTET+47) = 2
        NTBAR(3,NTET+47) = NTET+48
        NTSEITE(3,NTET+47) = 4
        NTBAR(4,NTET+47) = NTET+46
        NTSEITE(4,NTET+47) = 3
      ELSE
        NTBAR(1:4,NTET+47) = -1
        NTSEITE(1:4,NTET+47) = -1
      END IF


      NTECK(1,NTET+48) = INDCO(22)
      NTECK(2,NTET+48) = INDCO(19)
      NTECK(3,NTET+48) = INDCO(23)
      NTECK(4,NTET+48) = INDCO(14)
      CALL EINFUEGEN (INDCO(22),NTET+48)
      CALL EINFUEGEN (INDCO(19),NTET+48)
      CALL EINFUEGEN (INDCO(23),NTET+48)
      CALL EINFUEGEN (INDCO(14),NTET+48)
      RTCEN(NTET+48) = 0.25_DP* (FEM_KOOR(1,22)+FEM_KOOR(1,19)+
     .                           FEM_KOOR(1,23)+FEM_KOOR(1,14))
      STCEN(NTET+48) = 0.25_DP* (FEM_KOOR(2,22)+FEM_KOOR(2,19)+
     .                           FEM_KOOR(2,23)+FEM_KOOR(2,14))
      TTCEN(NTET+48) = 0.25_DP* (FEM_KOOR(3,22)+FEM_KOOR(3,19)+
     .                           FEM_KOOR(3,23)+FEM_KOOR(3,14))

      IF (COORD_TEST(INDCO(22),INDCO(19),INDCO(23),INDCO(14))) THEN
        NTBAR(2,NTET+48) = NTET+29
        NTSEITE(2,NTET+48) = 2
        NTBAR(3,NTET+48) = NTET+41
        NTSEITE(3,NTET+48) = 4
        NTBAR(4,NTET+48) = NTET+47
        NTSEITE(4,NTET+48) = 3
      ELSE
        NTBAR(1:4,NTET+48) = -1
        NTSEITE(1:4,NTET+48) = -1
      END IF


      DO ITET = NTET+1, NTET+48
        DO IS = 1,4
          JTET = NTBAR(IS,ITET)         ! NUMBER OF NEIGHBORING TETRAHEDRON
          IF (JTET > 0) THEN
            JS = NTSEITE(IS,ITET)
            IF (NTBAR(JS,JTET) < 0) THEN  ! SIDE POINTS TO A COLLAPSED
              NTBAR(IS,ITET) = 0          ! TETRAHEDRON
              NTSEITE(IS,ITET) = 0
            END IF
          END IF
        END DO
      END DO

      NTET = NTET+48

      RETURN

      CONTAINS

      SUBROUTINE EINFUEGEN (IC,ITET)
        INTEGER, INTENT(IN) :: IC, ITET
        TYPE(TET_ELEM), POINTER :: CUR

        ALLOCATE (CUR)
        CUR%NOTET = ITET
        CUR%NEXT_TET => COORTET(IC)%PTET
        COORTET(IC)%PTET => CUR
        MCLSTR = MCLSTR+1
      END SUBROUTINE EINFUEGEN


      FUNCTION COORD_TEST (I1,I2,I3,I4)
        INTEGER, INTENT(IN) :: I1,I2,I3,I4
        LOGICAL COORD_TEST
        LOGICAL LTEST

        LTEST= (I1==I2) .OR. (I1==I3) .OR. (I1==I4) .OR.
     .         (I2==I3) .OR. (I2==I4) .OR. (I3==I4)
        COORD_TEST = .NOT. LTEST
        RETURN
      END FUNCTION COORD_TEST

      END

C ===== SOURCE: suche_nachbarn.f


      SUBROUTINE SUCHE_NACHBARN

      USE CTETRA

      IMPLICIT NONE

      TYPE(TET_ELEM), POINTER :: CUR, CUR2
      INTEGER :: ITET,IS,JTET,JS,MINIS,MAXIS,MITIS,MINJS,MAXJS,MITJS,
     .           IP1,IP2,IP3,JP1,JP2,JP3, IC, i, j
      INTEGER :: JP(3),ip(3)
      INTEGER :: ITSIDE(3,4)
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/

      DO ITET=1,NTET      ! FOR ALL TETRAHEDRONS
        DO IS=1,4         ! AND FOR ALL SIDES OF EACH TETRAHEDRON
          IF (NTBAR(IS,ITET) == 0) THEN   ! IF IT HAS NO NEIGHBOR JET
            IP(1)=NTECK(ITSIDE(1,IS),ITET)
            IP(2)=NTECK(ITSIDE(2,IS),ITET)
            IP(3)=NTECK(ITSIDE(3,IS),ITET)

            CUR => COORTET(IP(1))%PTET
            WHLOOP:DO WHILE (ASSOCIATED(CUR))
              JTET = CUR%NOTET
              IF (JTET /= ITET) THEN  ! OMIT TETRAHEDRON ITET
                JSLOOP:DO JS=1,4                      ! CHECK ALL SIDES
                  IF (NTBAR(JS,JTET) == 0) THEN
                    JP(1)=NTECK(ITSIDE(1,JS),JTET)
                    JP(2)=NTECK(ITSIDE(2,JS),JTET)
                    JP(3)=NTECK(ITSIDE(3,JS),JTET)
                     iloop:do i=1,3
                      jloop:do j=i,3
                        if (ip(i) == jp(j)) then
                          ip1=jp(j)
                          jp(j)=jp(i)
                          jp(i)=ip1
                          cycle iloop
                        endif
                      end do jloop
                      cycle jsloop
                      end do iloop
                          NTBAR(IS,ITET) = JTET ! NEIGHBOR FOUND
                          NTSEITE(IS,ITET) = JS
                          NTBAR(JS,JTET) = ITET
                          NTSEITE(JS,JTET) = IS
                          EXIT WHLOOP
                  END IF
                END DO JSLOOP ! JS
              END IF
              CUR => CUR%NEXT_TET
            END DO WHLOOP ! WHILE
          END IF
        END DO
      END DO

!      DO IC=1,NCOOR
!        CUR => COORTET(IC)%PTET
!        DO WHILE (ASSOCIATED(CUR))
!           CUR2 => CUR
!           CUR => CUR%NEXT_TET
!           DEALLOCATE (CUR2)
!        END DO
!        NULLIFY(COORTET(IC)%PTET)
!      END DO
      WRITE (55,'(A,T25,I15)') ' Nachbar-Liste ',MCLSTR*8

      RETURN
      END







C ===== SOURCE: tet_step.f


      SUBROUTINE TET_STEP (IS,ITET,ISIDE,NRS)

      USE PRECISION
      USE COMUSR
      USE CTETRA
      USE CSTEP

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IS, ITET, ISIDE
      INTEGER, INTENT(INOUT) :: NRS
      REAL(DP) :: ARTRI3
      INTEGER :: I1, I2, I3
      INTEGER ITSIDE(3,4)
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/

! SIDE ISIDE OF TETRAHEDRON ITET BELONGS TO STEPFUNCTION ISTEP

      NRS=NRS+1
      I1=NTECK(ITSIDE(1,ISIDE),ITET)
      I2=NTECK(ITSIDE(2,ISIDE),ITET)
      I3=NTECK(ITSIDE(3,ISIDE),ITET)
      RRSTEP(IS,NRS+1)=RRSTEP(IS,NRS) +
     .               ARTRI3(XTETRA(I1),YTETRA(I1),ZTETRA(I1),
     .                      XTETRA(I2),YTETRA(I2),ZTETRA(I2),
     .                      XTETRA(I3),YTETRA(I3),ZTETRA(I3))
      IRSTEP(IS,NRS)=ITET
      IPSTEP(IS,NRS)=ISIDE
      ITSTEP(IS,NRS)=1
      IASTEP(IS,NRS)=0
      IBSTEP(IS,NRS)=1

      RETURN
      END
C ===== SOURCE: tetra_schnitt.f



      subroutine TETRA_schnitt
     .  (ebene, tetra, spanz, spunkt)

! Funktion zur Berechnung saemtlicher Schnittpunkte einer Ebene und
! eines Tetraeders (max.4).
! - ebene ist ein Vektor mit den Parametern (a,b,c,d) der Ebene
!   (ax+by+cz+d=0)
! - tetra ist eine 4x3-Matrix mit den (3-dim.) Eckpunkten des Tetraeders
! - spanz gibt die Anzahl der gefundenen Schnittpunkte an
! - in spunkt sind die gefundenen Schnittpunkte gespeichert

        USE PRECISION
        USE CCONA
        implicit none

        real(dp), intent(in), dimension(4)       :: ebene
        real(dp), intent(in), dimension(4,3)     :: tetra
        integer, intent(out)                     :: spanz
        real(dp), intent(out), dimension(4,3)    :: spunkt
        real(dp),dimension(3)                    :: punkt,SPP,teti,tetj
        real(dp),dimension(5,5)                  :: dmat
        REAL(DP)                                 :: DST, DSTMIN
        integer                                  :: i, j, k, l, ip
        INTEGER,dimension(2)                     :: ipdst
        logical                                  :: exi_spunkt

        spanz=0; i=1; j=1
        DSTMIN=HUGE(1.D0)
        dmat=DSTMIN
    !Schleife ueber alle moeglichen Kombinationen der Tetra-Eckpunkte
    !als Punkte einer Geraden
        do i=1,3
          j=(i+1)
          teti = tetra(i,:)
          do
            if (j==5) then
              exit
            end if
            tetj = tetra(j,:)
          !Berechnung eines Schnittpunktes der Ebene und einer
          !Kante des Tetraeders
!           call schnitt_ger_eb(spunkt(spanz+1,:), exi_spunkt,
!    .                          ebene, teti, tetj)
            call schnitt_ger_eb(spp, exi_spunkt,
     .                          ebene, teti, tetj)


          !Abfrage, ob Schnittpunkt gefunden
            if (exi_spunkt) then
               DO IP=1,SPANZ
                 dmat(ip,spanz+1)=SQRT(SUM((SPUNKT(IP,:)-SPP)**2))
               END DO
               dstmin=minval(dmat)
               IF (spanz < 4) then
                 if (minval(dmat(1:spanz,spanz+1)) > eps5) THEN
                   spanz = spanz+1
                   SPUNKT(SPANZ,:) = SPP
                 END IF
               ELSE
                 ipdst=minloc(dmat)
                 IF (IPDST(2) .NE. SPANZ+1) THEN
                   SPUNKT(IPDST(2),:) = SPP
                   dmat(1:ipdst(2)-1,ipdst(2))=
     .                  dmat(1:ipdst(2)-1,spanz+1)
                   dmat(ipdst(2),ipdst(2)+1:spanz)=
     .                  dmat(ipdst(2)+1:spanz,spanz+1)
                 END IF
               END IF
            end if

            j = j+1
          end do
        end do


    !bei vier Schnittpunkten
        if (spanz == 4) then

          call sort_ueberpruef(spunkt)

        end if

      end subroutine TETRA_schnitt



!----------------------------------------------------------------------

      subroutine schnitt_ger_eb(schnittpunkt,schn_ctrl,a,x0,x1)

! Funktion zur Berechnung eines Schnittpunktes von einer Ebene
! und einer Kante eines Tetraeders:
! Wenn ein Schnittpunkt vorhanden ist, ist schn_ctrl = .true.
! und der Schnittpunkt ist in schnittpunkt gespeichert.
! - a ist der Vektor mit den Koeffizienten der Ebene
! - x0 und x1 sind die Eckpunkte der Kanten des Tetraeders

        USE PRECISION
        USE CCONA
        implicit none

        logical, intent(out)               :: schn_ctrl
        real(dp), intent(in), dimension(4)   :: a
        real(dp), intent(in), dimension(3)   :: x0, x1
        real(dp), intent(out), dimension(3)   :: schnittpunkt
        real(dp),dimension(3)   :: x
        real(dp),dimension(3)   :: r, rnorm
        real(dp)                :: schnittfakt, hilfe
        real(dp)                :: betrag0, betrag1, betrag

        schn_ctrl = .true.

        !Berechnung des Richtungsvektors der Geraden,
        !die durch die Kante des Tetraeders beschrieben ist
        r = x1 - x0

        !Skalarprodukt von Normalenvektor der Ebene und
        !Richtungsvektor der Geraden
        hilfe = dot_product(a(2:4),r)

        !Abstand von x0 zur Ebene
        betrag0 = a(1) + dot_product(a(2:4),x0)

        !Abstand von x1 zur Ebene
        betrag1 = a(1) + dot_product(a(2:4),x1)

        !Abstand zwischen x0 und x1
        betrag = sqrt (sum (r * r))

        !normierter Richtungsvektor der Geraden
        rnorm = r / betrag

        schnittfakt = huge(1._dp)

        !Kante und Ebene sind parallel -> kein Schnittpunkt
!        if (abs(hilfe/(betraga*betragr)) < 1.E-4_DP) then
!        if (abs(hilfe) < 1.E-4_DP) then
!         if ((abs(hilfe) < eps5).or.
!     .       (betrag0*betrag1 > eps10) .or.
        if ((abs(dot_product(a(2:4),rnorm)) < eps5).or.
     .       (abs(betrag0-betrag1)/betrag < eps5)) then
           schn_ctrl = .false.
        else
           !Berechnung des Parameters der Geraden, der den
           !Schnittpunkt festlegt
           schnittfakt = -( (a(1) + dot_product(a(2:4),x0)) / hilfe )
           if (abs(schnittfakt) < eps5) schnittfakt=0._dp
        end if


        !Abfrage, ob der Schnittpunkt auch auf der Kante des Tetraeders liegt
        if ((schnittfakt >= 0._dp) .and.
     .      (schnittfakt <= 1._dp+eps5)) then
           !Berechnung des Schnittpunktes der Geraden und der Ebene
           x = x0 + schnittfakt * r
           schnittpunkt = x
        else
           schn_ctrl = .false.
        end if

      end subroutine schnitt_ger_eb


      subroutine sort_ueberpruef(punkt)

         USE PRECISION
         USE CCONA
         implicit none

         real(dp), intent(inout), dimension(4,3) :: punkt
         real(dp), dimension(3,2)                :: a
         real(dp), dimension(3)                  :: b, t, cen
         real(dp)                                :: x1, x2
         real(dp)                                :: d1, d2, d3, d4, d
         logical                                 :: lsw1, lsw2

         cen(1)=sum(punkt(1:4,1))*0.25_dp
         cen(2)=sum(punkt(1:4,2))*0.25_dp
         cen(3)=sum(punkt(1:4,3))*0.25_dp

         d1 = sqrt(sum((punkt(1,:)-cen)**2))
         d2 = sqrt(sum((punkt(2,:)-cen)**2))
         d3 = sqrt(sum((punkt(3,:)-cen)**2))
         d4 = sqrt(sum((punkt(4,:)-cen)**2))
         d = 1._dp/max(min(d1,d2,d3,d4),eps10)

         lsw1=.true.
         lsw2=.true.
         do while (lsw1 .or. lsw2)

!  test p1-p2 and p3-p4
            a(1,1) = punkt(2,1) - punkt(1,1)
            a(2,1) = punkt(2,2) - punkt(1,2)
            a(3,1) = punkt(2,3) - punkt(1,3)
            a(1,2) = punkt(3,1) - punkt(4,1)
            a(2,2) = punkt(3,2) - punkt(4,2)
            a(3,2) = punkt(3,3) - punkt(4,3)

            b(1) = punkt(3,1) - punkt(1,1)
            b(2) = punkt(3,2) - punkt(1,2)
            b(3) = punkt(3,3) - punkt(1,3)

            a = a*d
            b = b*d

            call solve_lgs

            if ((x1 >= 0._dp) .and. (x1 <= 1._dp) .and.
     .          (x2 >= 0._dp) .and. (x2 <= 1._dp)) then
!  intersection found, switch points 2 and 3
               t(1:3) = punkt(2,1:3)
               punkt(2,1:3) = punkt(3,1:3)
               punkt(3,1:3) = t(1:3)
               lsw1 = .true.
            else
               lsw1 = .false.
            end if

!  test p1-p4 and p2-p3
            a(1,1) = punkt(4,1) - punkt(1,1)
            a(2,1) = punkt(4,2) - punkt(1,2)
            a(3,1) = punkt(4,3) - punkt(1,3)
            a(1,2) = punkt(2,1) - punkt(3,1)
            a(2,2) = punkt(2,2) - punkt(3,2)
            a(3,2) = punkt(2,3) - punkt(3,3)

            b(1) = punkt(2,1) - punkt(1,1)
            b(2) = punkt(2,2) - punkt(1,2)
            b(3) = punkt(2,3) - punkt(1,3)

            a = a*d
            b = b*d

            call solve_lgs

            if ((x1 >= 0._dp) .and. (x1 <= 1._dp) .and.
     .          (x2 >= 0._dp) .and. (x2 <= 1._dp)) then
!  intersection found, switch points 3 and 4
               t(1:3) = punkt(3,1:3)
               punkt(3,1:3) = punkt(4,1:3)
               punkt(4,1:3) = t(1:3)
               lsw2 = .true.
            else
               lsw2 = .false.
            end if

         end do                 ! do while

         contains

         subroutine solve_lgs
           real(dp) :: c(2,2), d(2)
           real(dp) :: detc, detc1, detc2

           c(1,1) = a(1,1)**2 + a(2,1)**2 + a(3,1)**2
           c(1,2) = a(1,1)*a(1,2) + a(2,1)*a(2,2) + a(3,1)*a(3,2)
           c(2,1) = c(1,2)
           c(2,2) = a(1,2)**2 + a(2,2)**2 + a(3,2)**2

           d(1) = a(1,1)*b(1) + a(2,1)*b(2) + a(3,1)*b(3)
           d(2) = a(1,2)*b(1) + a(2,2)*b(2) + a(3,2)*b(3)

           detc  = c(1,1)*c(2,2) - c(1,2)*c(2,1)
           detc1 = d(1)*c(2,2) - d(2)*c(1,2)
           detc2 = c(1,1)*d(2) - c(2,1)*d(1)

           if (abs(detc) < 1.e-6_dp) then
             x1 = 10._dp
             x2 = 10._dp
           else
             x1 = detc1/detc
             x2 = detc2/detc
           end if

       end subroutine solve_lgs

       end subroutine sort_ueberpruef


