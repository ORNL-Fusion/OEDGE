

      INTEGER FUNCTION LEARCT (X,Y,Z)

      USE PRECISION
      USE COMUSR
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

      WRITE (6,*) ' POINT ',X,Y,Z
      WRITE (6,*) ' OUTSIDE OF ALL TETRAHEDRONS '
      LEARCT=0
      RETURN
      END


      INTEGER FUNCTION LEARCT_old (X,Y,Z)

      USE PRECISION
      USE COMUSR
      USE CCONA
      USE CTETRA

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X,Y,Z
      REAL(DP) :: PC1(3), PC2(3), PC3(3), PC4(3), P(3)
      REAL(DP) :: V1, V2, V3, V4, CAL_VOL
      INTEGER :: ITET

      P(1:3) = (/ X, Y, Z /)

      DO ITET=1,NTET
        IF (VOL(ITET) > 1.E-6) THEN
          PC1(1:3)= (/ XTETRA(NTECK(1,ITET)), YTETRA(NTECK(1,ITET)),
     .                 ZTETRA(NTECK(1,ITET)) /)
          PC2(1:3)= (/ XTETRA(NTECK(2,ITET)), YTETRA(NTECK(2,ITET)),
     .                 ZTETRA(NTECK(2,ITET)) /)
          PC3(1:3)= (/ XTETRA(NTECK(3,ITET)), YTETRA(NTECK(3,ITET)),
     .                 ZTETRA(NTECK(3,ITET)) /)
          PC4(1:3)= (/ XTETRA(NTECK(4,ITET)), YTETRA(NTECK(4,ITET)),
     .                 ZTETRA(NTECK(4,ITET)) /)

          IF ((MIN(PC1(1),PC2(1),PC3(1),PC4(1)) <= X) .AND.
     .        (MAX(PC1(1),PC2(1),PC3(1),PC4(1)) >= X) .AND.
     .        (MIN(PC1(2),PC2(2),PC3(2),PC4(2)) <= Y) .AND.
     .        (MAX(PC1(2),PC2(2),PC3(2),PC4(2)) >= Y) .AND.
     .        (MIN(PC1(3),PC2(3),PC3(3),PC4(3)) <= Z) .AND.
     .        (MAX(PC1(3),PC2(3),PC3(3),PC4(3)) >= Z)) THEN

            V1 = CAL_VOL (PC1,PC2,PC3,P)
            V2 = CAL_VOL (PC3,PC2,PC4,P)
            V3 = CAL_VOL (PC1,PC3,PC4,P)
            V4 = CAL_VOL (PC1,PC4,PC2,P)

            IF ((ABS(V1+V2+V3+V4-VOL(ITET)) < 1.D-3*VOL(ITET)) .AND.
     .          (MIN(V1,V2,V3,V4) >= -EPS5*VOL(ITET))) THEN
              LEARCT_old=ITET
              RETURN
            END IF
          END IF
        END IF
      END DO

      WRITE (6,*) ' POINT ',X,Y,Z
      WRITE (6,*) ' OUTSIDE OF ALL TETRAHEDRONS '
      LEARCT_old=0
      RETURN
      END
