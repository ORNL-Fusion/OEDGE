      FUNCTION LEARC1 (X,Y,Z,IPO,IAN,IEN,LOGX,LOGY,NP,TEXT)
C
C  REV. JULY 01: LEVGEO=3: ALL CALCULATIONS IN RELATIVE DISTANCES
C  REV. JULY 01: LEVGEO=4: ALL CALCULATIONS IN RELATIVE DISTANCES
C
C
C   LOGX=TRUE: PARTICLE IS ON A RADIAL SURFACE
C   LOGY=TRUE: PARTICLE IS ON A POLOIDAL SURFACE
C
C   FIND RADIAL MESHPOINT NUMBER LEARC1,
C   (AND POLYGON INDEX IPO, IF NLPLG)
C  LEVGEO=3:
C   IF .NOT.LOGX AND .NOT.LOGY
C     SEARCH IN RADIAL CELLS IR: [IAN,IEN], I.E.
C     SEARCH BETWEEN (!!!) RADIAL SURFACES IAN AND IEN+1
C     THIS SEARCH COVERS THE WHOLE POLOIDAL RANGE
C   IF LOGX
C     SEARCH ON (!!!) RADIAL SURF. IAN FOR POLOIDAL MESH NUMBER IPO
C   IF LOGY
C     SEARCH ON (!!!) POLOIDAL SURF. IAN FOR RADIAL MESH NUMBER LEARC1
C  LEVGEO=4:
C
      USE PRECISION
      USE PARMMOD
      USE CGRID
      USE CGEOM
      USE CCONA
      USE CPOLYG
      USE CLOGAU
      USE CTRIG

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X, Y, Z
      INTEGER, INTENT(IN) :: IAN, IEN, NP
      INTEGER, INTENT(OUT) :: IPO
      CHARACTER(LEN=*), INTENT(IN) :: TEXT
      LOGICAL, INTENT(IN) ::  LOGX, LOGY

      REAL(DP) :: X1, X2, X3, X4, Y1, Y2, Y3, Y4, DET4, XCMIN, XCMAX,
     .          ERRMIN, YCMIN, YCMAX, ERR7, YQ, DX4, ERR4, ATQ, XE, XEQ,
     .          XMX1, XMX2, YMY1, YMY2, YMY4, ERR1, DX1, XMX4, DET3,
     .          DET1, DET2, D1, D2, D3, XTRMIN, XTRMAX, YTRMIN, YTRMAX, 
     .          D4, DELTAX, DELTAY, O23, O12, O31, O41, O34, O13, XS1,
     .          XS2, YS1, YS2, O1, O2
      REAL(DP), SAVE :: XMIN, YMIN, DISTX, DISTY, XMAX, YMAX, 
     .                  EPDY, EPDXDY, EPDX
      INTEGER, SAVE :: IFIRST
      INTEGER :: K, L, IM, LM, IEP, KH, LEARCT, LEAUSR, IMARK, LMARK,
     .           I, J, IE, LEARC1, IA, INTR1, INTR2,
     .           IX, IY, INUM, IHEADX1, IHEADX2, IHEADY1, IHEADY2

      REAL(DP), ALLOCATABLE, SAVE ::
     R D12(:,:), D12I(:,:), D14(:,:), D14I(:,:), OBSC(:,:)
      LOGICAL :: LG(N1ST,N2NDPLG), LG1(NRAD)

CTK DATENSTRUKTUR FUER DREIECKS UND VIERECKSGITTER
      TYPE :: CELL
        INTEGER :: TRIANGLE
        TYPE(CELL),POINTER :: NEXT
      END TYPE CELL

      TYPE :: CELL4
        INTEGER :: IX
        INTEGER :: IY
        TYPE(CELL4),POINTER :: NEXT
      END TYPE CELL4

      TYPE :: POIFELD
        TYPE (CELL),POINTER :: P
      END TYPE POIFELD

      TYPE :: POI4
        TYPE (CELL4),POINTER :: P
      END TYPE POI4

      TYPE (POIFELD) :: HELPCUR(4)
      TYPE (POIFELD),ALLOCATABLE,SAVE :: HEADS(:,:)
      TYPE (CELL),POINTER :: CUR
      TYPE (POI4) :: HELPCUR4(4)
      TYPE (POI4),ALLOCATABLE,SAVE :: HEADS4(:,:)
      TYPE (CELL4),POINTER :: CUR4,HELPP
C
       
      DATA IFIRST /0/
C
      IA=IAN
      IE=IEN
C
      IF (LEVGEO.EQ.4) THEN
C
        IF (IFIRST.EQ.0) THEN
          ALLOCATE (OBSC(N1ST,1))
          IFIRST = 1
          ALLOCATE(HEADS(100,100))
          DO I=1,100
            DO J=1,100
              NULLIFY(HEADS(I,J)%P)
            ENDDO
          ENDDO
          XMIN=MINVAL(XTRIAN(1:NRKNOT))
          YMIN=MINVAL(YTRIAN(1:NRKNOT))
          XMAX=MAXVAL(XTRIAN(1:NRKNOT))
          YMAX=MAXVAL(YTRIAN(1:NRKNOT))
          XMIN = MIN( XMIN *(1.-EPS5), XMIN *(1.+EPS5) )
          XMAX = MAX( XMAX *(1.-EPS5), XMAX *(1.+EPS5) )
          YMIN = MIN( YMIN *(1.-EPS5), YMIN *(1.+EPS5) )
          YMAX = MAX( YMAX *(1.-EPS5), YMAX *(1.+EPS5) )
          DISTX=(XMAX-XMIN)/100._DP
          DISTY=(YMAX-YMIN)/100._DP
          EPDX=DISTX*EPS5
          EPDY=DISTY*EPS5
          EPDXDY=EPDX*EPDY
          DO I=1,NTRII
C  CHECK CELLS FOR ORIENTATION, DISTORTION, CONVEX SHAPE, ETC...
              X1=XTRIAN(NECKE(1,I))
              Y1=YTRIAN(NECKE(1,I))
              X2=XTRIAN(NECKE(2,I))
              Y2=YTRIAN(NECKE(2,I))
              X3=XTRIAN(NECKE(3,I))
              Y3=YTRIAN(NECKE(3,I))
C  CENTER OF GRAVITY IN TRIANGLE 1 AND TRIANGLE 2
              XS1=(X1+X2+X3)/3._DP
              YS1=(Y1+Y2+Y3)/3._DP
              D1 =(X1-XS1)*(Y2-YS1)-(Y1-YS1)*(X2-XS1)
              D2 =(X2-XS1)*(Y3-YS1)-(Y2-YS1)*(X3-XS1)
              D3 =(X3-XS1)*(Y1-YS1)-(Y3-YS1)*(X1-XS1)
              DET1=MIN(0._DP,D1+EPDXDY)
              DET2=MIN(0._DP,D2+EPDXDY)
              DET3=MIN(0._DP,D3+EPDXDY)
              O1=-1.
              IF (ABS(DET1+DET2+DET3) .GE. EPDXDY) THEN
C  SINCE WE KNOW THAT S1 IS INSIDE, THIS
C  MUST BE DUE TO OPPOSITE ORIENTATION OF THIS TRIANGLE
                O1=1.
              ENDIF
C
              IF (O1.LT.0) THEN
C  NORMAL CASE: ORIENTATION IN TRIANGLE IS NEGATIVE (MATHEM. POSITIVE)
                OBSC(I,1)=0._DP
              ELSEIF (O1.GT.0) THEN
C  OPPOSITE GRID ORIENTATION: ORIENTATION IN TRIANGLE IS POSITIVE
                OBSC(I,1)=-1._DP
                WRITE (6,*) 'OBSCURE CELL DETECTED, ITRII= ',I
                WRITE (6,*) 'OPPOSITE ORIENTATION IN TRIANGLE '
              ENDIF
          ENDDO
C
          DO I=1,NTRII
            XTRMIN = MIN(XTRIAN(NECKE(1,I)),XTRIAN(NECKE(2,I)),
     .                   XTRIAN(NECKE(3,I)))
            XTRMAX = MAX(XTRIAN(NECKE(1,I)),XTRIAN(NECKE(2,I)),
     .                   XTRIAN(NECKE(3,I)))
            YTRMIN = MIN(YTRIAN(NECKE(1,I)),YTRIAN(NECKE(2,I)),
     .                   YTRIAN(NECKE(3,I)))
            YTRMAX = MAX(YTRIAN(NECKE(1,I)),YTRIAN(NECKE(2,I)),
     .                   YTRIAN(NECKE(3,I)))
            DELTAX=XTRMIN-XMIN
            IHEADX1=INT(DELTAX/DISTX)+1
            DELTAX=XTRMAX-XMIN
            IHEADX2=MAX(IHEADX1,INT(DELTAX/DISTX)+1)
            DELTAY=YTRMIN-YMIN
            IHEADY1=INT(DELTAY/DISTY)+1
            DELTAY=YTRMAX-YMIN
            IHEADY2=MAX(IHEADY1,INT(DELTAY/DISTY)+1)
            DO IX=IHEADX1,IHEADX2
              DO IY=IHEADY1,IHEADY2
                ALLOCATE(CUR)
                CUR%TRIANGLE = I
                CUR%NEXT => HEADS(IX,IY)%P
                HEADS(IX,IY)%P => CUR
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
C  END OF IFIRST SEGMENT FOR TRIANGELS
C
        INUM=0
C
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

        HELPCUR(1)%P => HEADS(IHEADX1,IHEADY1)%P
        IF (IHEADX2 .GT. 0) THEN
          HELPCUR(2)%P => HEADS(IHEADX2,IHEADY1)%P
        ELSE
          NULLIFY(HELPCUR(2)%P)
        ENDIF
        IF (IHEADY2 .GT. 0) THEN
          HELPCUR(3)%P => HEADS(IHEADX1,IHEADY2)%P
        ELSE
          NULLIFY(HELPCUR(3)%P)
        ENDIF
        IF ((IHEADX2 .GT. 0) .AND. (IHEADY2 .GT. 0)) THEN
          HELPCUR(4)%P => HEADS(IHEADX2,IHEADY2)%P
        ELSE
          NULLIFY(HELPCUR(4)%P)
        ENDIF

        LG1=.FALSE.
        DO J=1,4
          DO WHILE (ASSOCIATED(HELPCUR(J)%P))
            I = HELPCUR(J)%P%TRIANGLE
C  CELL I ALREADY TESTED BEFORE ?
            IF (LG1(I)) GOTO 5
            LG1(I)=.TRUE.
            D1 = (XTRIAN(NECKE(1,I))-X)*(YTRIAN(NECKE(2,I))-Y)-
     .           (YTRIAN(NECKE(1,I))-Y)*(XTRIAN(NECKE(2,I))-X)
            D2 = (XTRIAN(NECKE(2,I))-X)*(YTRIAN(NECKE(3,I))-Y)-
     .           (YTRIAN(NECKE(2,I))-Y)*(XTRIAN(NECKE(3,I))-X)
            D3 = (XTRIAN(NECKE(3,I))-X)*(YTRIAN(NECKE(1,I))-Y)-
     .           (YTRIAN(NECKE(3,I))-Y)*(XTRIAN(NECKE(1,I))-X)
            DET1=MIN(0._DP,D1+EPDXDY)
            DET2=MIN(0._DP,D2+EPDXDY)
            DET3=MIN(0._DP,D3+EPDXDY)
            IF (ABS(DET1+DET2+DET3) .LT. EPDXDY) THEN
              INUM=INUM+1
              IM = I
            ENDIF
5           HELPCUR(J)%P => HELPCUR(J)%P%NEXT
          ENDDO
        ENDDO
        IF (IM.LT.1.OR.IM.GT.NTRII) THEN
          WRITE (6,*) 'NO TRIANGLE FOUND IN LEARC1 FOR '
          WRITE (6,*) 'X = ',X,' Y = ',Y
          WRITE (6,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
          WRITE (6,*) 'NPANU,IM= ',NP,IM
        ELSEIF (INUM.GT.1) THEN
          CALL MASAGE ('WARNING FROM LEARC1, INUM.GT.1               ')
          CALL MASR2('X,Y             ',X,Y)
          WRITE (6,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
          WRITE (6,*) 'NPANU,INUM,IM= ',NP,INUM,IM
          WRITE (6,*) 'IAN,IEN,LOGX,LOGY ',IAN,IEN,LOGX,LOGY
        ENDIF
        LEARC1=IM
C
      ELSEIF (LEVGEO.EQ.3) THEN
C
        IF (IFIRST.EQ.0) THEN
          IFIRST=1
          ALLOCATE (D12(N1ST,N2ND))
          ALLOCATE (D12I(N1ST,N2ND))
          ALLOCATE (D14(N1ST,N2ND))
          ALLOCATE (D14I(N1ST,N2ND))
          ALLOCATE (OBSC(N1ST,N2ND))
          DO 1 I=1,NR1ST
            DO 2 L=1,NP2NDM
              D12(I,L)=SQRT((XPOL(I,L)-XPOL(I,L+1))**2+
     .                      (YPOL(I,L)-YPOL(I,L+1))**2)
              D12I(I,L)=1./(ABS(D12(I,L))+EPS60)
2           CONTINUE
1         CONTINUE
          DO 3 I=1,NR1STM
            DO 4 L=1,NP2ND
              D14(I,L)=SQRT((XPOL(I,L)-XPOL(I+1,L))**2+
     .                      (YPOL(I,L)-YPOL(I+1,L))**2)
              D14I(I,L)=1./(ABS(D14(I,L))+EPS60)
4           CONTINUE
3         CONTINUE
C
          ALLOCATE(HEADS4(100,100))
C  SET EQUIDISTANT X-Y GRID, WHICH COVERS POLYGON GRID
          DO IX=1,100
            DO IY=1,100
              NULLIFY(HEADS4(IX,IY)%P)
            ENDDO
          ENDDO
          XMIN=1.D60
          YMIN=1.D60
          XMAX=-1.D60
          YMAX=-1.D60
          DO I=1,NR1ST
            DO L=1,NP2ND
              XMIN = MIN(XMIN,XPOL(I,L))
              YMIN = MIN(YMIN,YPOL(I,L))
              XMAX = MAX(XMAX,XPOL(I,L))
              YMAX = MAX(YMAX,YPOL(I,L))
            ENDDO
          ENDDO
          XMIN = MIN( XMIN *(1._DP-EPS5), XMIN *(1._DP+EPS5) )
          XMAX = MAX( XMAX *(1._DP-EPS5), XMAX *(1._DP+EPS5) )
          YMIN = MIN( YMIN *(1._DP-EPS5), YMIN *(1._DP+EPS5) )
          YMAX = MAX( YMAX *(1._DP-EPS5), YMAX *(1._DP+EPS5) )
          DISTX=(XMAX-XMIN)/100.
          DISTY=(YMAX-YMIN)/100.
          EPDX=DISTX*EPS10
          EPDY=DISTY*EPS10
          EPDXDY=EPDX+EPDY
C  FOR EACH POLYGON CELL (I,L) FIND THE RANGE IHEADX1,....IHEADY2
C                              SUCH THAT THIS CELL (I,L) IS ENTIRELY
C                              IN THAT SECTION OF THE REGULAR IX,IY GRID
          DO I=1,NR1STM
            DO K=1,NPPLG
            DO L=NPOINT(1,K),NPOINT(2,K)-1
              XCMIN=MIN(XPOL(I,L),XPOL(I+1,L),XPOL(I+1,L+1),XPOL(I,L+1))
              XCMAX=MAX(XPOL(I,L),XPOL(I+1,L),XPOL(I+1,L+1),XPOL(I,L+1))
              YCMIN=MIN(YPOL(I,L),YPOL(I+1,L),YPOL(I+1,L+1),YPOL(I,L+1))
              YCMAX=MAX(YPOL(I,L),YPOL(I+1,L),YPOL(I+1,L+1),YPOL(I,L+1))
              DELTAX=XCMIN-XMIN
              IHEADX1=INT(DELTAX/DISTX)+1
              DELTAX=XCMAX-XMIN
              IHEADX2=MAX(IHEADX1,INT(DELTAX/DISTX)+1)
              DELTAY=YCMIN-YMIN
              IHEADY1=INT(DELTAY/DISTY)+1
              DELTAY=YCMAX-YMIN
              IHEADY2=MAX(IHEADY1,INT(DELTAY/DISTY)+1)
C  ADD POLYGON CELL (I,L) TO THE LIST FOR ALL REGULAR (IX,IY) CELLS
C  IN THE RANGE IHEADX1,....IHEADY2
              DO IX=IHEADX1,IHEADX2
                DO IY=IHEADY1,IHEADY2
                  ALLOCATE(CUR4)
                  CUR4%IX = I
                  CUR4%IY = L
                  CUR4%NEXT => HEADS4(IX,IY)%P
                  HEADS4(IX,IY)%P => CUR4
                ENDDO
              ENDDO
C  CHECK CELLS FOR ORIENTATION, DISTORTION, CONVEX SHAPE, ETC...
              X1=XPOL(I,L)
              Y1=YPOL(I,L)
              X2=XPOL(I,L+1)
              Y2=YPOL(I,L+1)
              X3=XPOL(I+1,L+1)
              Y3=YPOL(I+1,L+1)
              X4=XPOL(I+1,L)
              Y4=YPOL(I+1,L)
C  FROM THIS POINT: QUADRANGLE IS SUBDIVIDED INTO 2 TRIANGLES
C  ALONG LINE 1 -- 3. THIS IS USED LATER IN LEARC1.
C
C  CENTER OF GRAVITY IN TRIANGLE 1 AND TRIANGLE 2
              XS1=(X1+X2+X3)/3._DP
              YS1=(Y1+Y2+Y3)/3._DP
              XS2=(X1+X4+X3)/3._DP
              YS2=(Y1+Y4+Y3)/3._DP
              D1 =(X1-XS1)*(Y2-YS1)-(Y1-YS1)*(X2-XS1)
              D2 =(X2-XS1)*(Y3-YS1)-(Y2-YS1)*(X3-XS1)
              D3 =(X3-XS1)*(Y1-YS1)-(Y3-YS1)*(X1-XS1)
              DET1=MIN(0._DP,D1+EPDXDY)
              DET2=MIN(0._DP,D2+EPDXDY)
              DET3=MIN(0._DP,D3+EPDXDY)
              O1=-1.
              IF (ABS(DET1+DET2+DET3) .GE. EPDXDY) THEN
C  SINCE WE KNOW THAT S1 IS INSIDE, THIS
C  MUST BE DUE TO OPPOSITE ORIENTATION OF THIS TRIANGLE
                O1=1.
              ENDIF
              D1 =-(X1-XS2)*(Y4-YS2)+(Y1-YS2)*(X4-XS2)
              D4 =-(X4-XS2)*(Y3-YS2)+(Y4-YS2)*(X3-XS2)
              D3 =-(X3-XS2)*(Y1-YS2)+(Y3-YS2)*(X1-XS2)
              DET1=MIN(0._DP,D1+EPDXDY)
              DET4=MIN(0._DP,D4+EPDXDY)
              DET3=MIN(0._DP,D3+EPDXDY)
              O2=-1.
              IF (ABS(DET1+DET4+DET3) .GE. EPDXDY) THEN
                O2=1.
              ENDIF
C
              IF (O1.LT.0.AND.O2.LT.0) THEN
C  NORMAL CASE: ORIENTATION IN BOTH TRIANGES IS NEGATIVE
                OBSC(I,L)=0._DP
              ELSEIF (O1.GT.0.AND.O2.GT.0) THEN
C  OPPOSITE GRID ORIENTATION: ORIENTATION IN BOTH TRIANGES IS POSITIVE
                OBSC(I,L)=-1._DP
              ELSE
C  OBSCURE CELL: CONVEX QUADRANGE, POINT 2 OR POINT 4 INSIDE
                WRITE (6,*) 'OBSCURE CELL DETECTED, IR,IP= ',I,L
                WRITE (6,*) 'PLASMA FIELD LIKELY TO BE CORRUPTED HERE '
                OBSC(I,L)=1._DP
              ENDIF
C  TO IDENTIFY OBSCURE CELLS, IN WHICH POINT 1 OR POINT 3 IS INSIDE,
C  THE SAME MUST BE REPEATED FOR THE QUADRANGE SUBDIVIDED ALONG
C  LINE 2--4 INTO TWO TRIANGLES.
C  FOR PURPOSES OF LEARC1 THESE CELLS NEED NOT BE IDENTIFIED, HOWEVER.
            ENDDO
            ENDDO
          ENDDO
        ENDIF
C
C  END OF IFIRST SEGMENT FOR QUADRANGELS
C
        INUM=0
C
        DELTAX=X-XMIN
        IHEADX2 = 0
        IF (ABS(MOD(DELTAX,DISTX)) .LT. EPDX) THEN
          IHEADX2=INT(DELTAX/DISTX)
        ENDIF
        IHEADX1=INT(DELTAX/DISTX)+1
C
        DELTAY=Y-YMIN
        IHEADY2 = 0
        IF (ABS(MOD(DELTAY,DISTY)) .LT. EPDY) THEN
          IHEADY2=INT(DELTAY/DISTY)
        ENDIF
        IHEADY1=INT(DELTAY/DISTY)+1
C
C  EITHER ONE, TWO, THREE OR FOUR REGULAR CELLS ARE CONSIDERED
C  THIS DEPENDS UPON WHERE THE POINT X,Y, LIES, RELATIVE TO THE
C  REGULAR GRID LINES
C
C  X,Y IS IN REGULAR CELL IX,IY=IHEADX1,IHEADY1. ALLWAYS CHECKED.
        HELPCUR4(1)%P => HEADS4(IHEADX1,IHEADY1)%P

C  X,Y POSSIBLY ALSO IN REGULAR CELL IHEADX2,IHEADY1?
        IF (IHEADX2 .GT. 0) THEN
          HELPCUR4(2)%P => HEADS4(IHEADX2,IHEADY1)%P
        ELSE
          NULLIFY(HELPCUR4(2)%P)
        ENDIF
C
C  X,Y POSSIBLY ALSO IN REGULAR CELL IHEADX1,IHEADY2?
        IF (IHEADY2 .GT. 0) THEN
          HELPCUR4(3)%P => HEADS4(IHEADX1,IHEADY2)%P
        ELSE
          NULLIFY(HELPCUR4(3)%P)
        ENDIF
C
C  X,Y POSSIBLY ALSO IN REGULAR CELL IHEADX2,IHEADY2?
        IF ((IHEADX2 .GT. 0) .AND. (IHEADY2 .GT. 0)) THEN
          HELPCUR4(4)%P => HEADS4(IHEADX2,IHEADY2)%P
        ELSE
          NULLIFY(HELPCUR4(4)%P)
        ENDIF
C
        ERRMIN=1.D30
        IF (LOGX) THEN
          IEP=IA
          GOTO 500
        ENDIF
        IF (LOGY) THEN
          IEP=IA
          GOTO 750
        ENDIF
C
        LG=.FALSE.
        DO J=1,4
          HELPP => HELPCUR4(J)%P
          DO WHILE (ASSOCIATED(HELPP))
            I = HELPP%IX
            L = HELPP%IY
C  CELL I,L ALREADY TESTED BEFORE ?
            IF (LG(I,L)) GOTO 20
            LG(I,L)=.TRUE.
            IF ((I.LT.IA) .OR. (I.GT.IE) .OR. (L.GT.NP2NDM)) GOTO 20
C  NORMAL CASE:  OBSC(I,L)=0
            IF (OBSC(I,L).EQ.0._DP) THEN
            X1=XPOL(I,L)
            Y1=YPOL(I,L)
            X2=XPOL(I,L+1)
            Y2=YPOL(I,L+1)
            X3=XPOL(I+1,L+1)
            Y3=YPOL(I+1,L+1)
C
            D1 =(X1-X)*(Y2-Y)-(Y1-Y)*(X2-X)
            D2 =(X2-X)*(Y3-Y)-(Y2-Y)*(X3-X)
            D3 =(X3-X)*(Y1-Y)-(Y3-Y)*(X1-X)
            DET1=MIN(0._DP,D1+EPDXDY)
            DET2=MIN(0._DP,D2+EPDXDY)
            DET3=MIN(0._DP,D3+EPDXDY)
            IF (ABS(DET1+DET2+DET3) .LT. EPDXDY) THEN
              INUM=INUM+1
              IM=I
              LM=L
              GOTO 20
            ENDIF
            X4=XPOL(I+1,L)
            Y4=YPOL(I+1,L)
            D1 =-(X1-X)*(Y4-Y)+(Y1-Y)*(X4-X)
            D4 =-(X4-X)*(Y3-Y)+(Y4-Y)*(X3-X)
            D3 =-D3
            DET1=MIN(0._DP,D1+EPDXDY)
            DET4=MIN(0._DP,D4+EPDXDY)
            DET3=MIN(0._DP,D3+EPDXDY)
            IF (ABS(DET1+DET4+DET3) .LT. EPDXDY) THEN
              INUM=INUM+1
              IM=I
              LM=L
              GOTO 20
            ENDIF
C  OPPOSITE ORIENTATION OF GRID:  OBSC(I,L)=-1
            ELSEIF (OBSC(I,L).EQ.-1._DP) THEN
              WRITE (6,*) 'OPPOSITE GRID ORIENTATION '
              WRITE (6,*) 'FUNCTION LEARC1 NEEDS TO BE EXTENDED'
              CALL EXIT_OWN(1)
C  OBSCURE CELL
C  THIS PART SHOULD ALREADY WORK FOR BOTH ORIENTATIONS OF THE GRID
            ELSEIF (OBSC(I,L).EQ.1._DP) THEN
            X1=XPOL(I,L)
            Y1=YPOL(I,L)
            X2=XPOL(I,L+1)
            Y2=YPOL(I,L+1)
            X3=XPOL(I+1,L+1)
            Y3=YPOL(I+1,L+1)
            X4=XPOL(I+1,L)
            Y4=YPOL(I+1,L)
C POINT MUST BE IN ONE TRIANGLE, BUT NOT IN THE OTHER AS WELL
            INTR1=0
            INTR2=0
            O12 =(X1-X)*(Y2-Y)-(Y1-Y)*(X2-X)
            O23 =(X2-X)*(Y3-Y)-(Y2-Y)*(X3-X)
            O31 =(X3-X)*(Y1-Y)-(Y3-Y)*(X1-X)
C  NORMAL ORIENTATION: ALL D POSITIVE ?
            DET1=MIN(0._DP,O12+EPDXDY)
            DET2=MIN(0._DP,O23+EPDXDY)
            DET3=MIN(0._DP,O31+EPDXDY)
            IF (ABS(DET1+DET2+DET3) .LT. EPDXDY) THEN
              INTR1=1
              IM=I
              LM=L
              GOTO 25
            ENDIF
C  OPPOSITE ORIENTATION: ALL D NEGATIVE ?
            DET1=MAX(0._DP,O12+EPDXDY)
            DET2=MAX(0._DP,O23+EPDXDY)
            DET3=MAX(0._DP,O31+EPDXDY)
            IF (ABS(DET1+DET2+DET3) .LT. EPDXDY) THEN
              INTR1=1
              IM=I
              LM=L
              GOTO 25
            ENDIF
25          CONTINUE
C   CHECK SECOND TRIANGLE
            O13 =-O31
            O34 =-(X4-X)*(Y3-Y)+(Y4-Y)*(X3-X)
            O41 =-(X1-X)*(Y4-Y)+(Y1-Y)*(X4-X)
C  NORMAL ORIENTATION: ALL D POSITIVE ?
            DET1=MIN(0._DP,O13+EPDXDY)
            DET4=MIN(0._DP,O34+EPDXDY)
            DET3=MIN(0._DP,O41+EPDXDY)
            IF (ABS(DET1+DET4+DET3) .LT. EPDXDY) THEN
              INTR2=1
              IM=I
              LM=L
              GOTO 27
            ENDIF
C  OPPOSITE ORIENTATION: ALL D NEGATIVE ?
            DET1=MAX(0._DP,O13+EPDXDY)
            DET4=MAX(0._DP,O34+EPDXDY)
            DET3=MAX(0._DP,O41+EPDXDY)
            IF (ABS(DET1+DET4+DET3) .LT. EPDXDY) THEN
              INTR2=1
              IM=I
              LM=L
              GOTO 27
            ENDIF
27          CONTINUE
            IF (INTR1+INTR2.EQ.1) THEN
              INUM=INUM+1
              GOTO 20
            ENDIF
            ENDIF
C
20          HELPP => HELPP%NEXT
          ENDDO
        ENDDO
C
        IF (INUM.EQ.1) GOTO 1000
        IEP=IE+1
C
C  CHECK FOR NEAREST BOUNDARY, BECAUSE NO VALID CELL INDEX FOUND
C  FIRST TRY RADIAL SURFACES
C  THIS SECTION ALSO: IF LOGX, CHECK ON RADIAL SURFACE IA
500     CONTINUE
        DO J=1,4
          HELPP => HELPCUR4(J)%P
          DO WHILE (ASSOCIATED(HELPP))
            I = HELPP%IX
            L = HELPP%IY
121         IF (L.GT.NP2NDM) GOTO 21
            XMX1=X-XPOL(I,L)
            YMY1=Y-YPOL(I,L)
            XMX2=X-XPOL(I,L+1)
            YMY2=Y-YPOL(I,L+1)
            D1=SQRT(XMX1*XMX1+YMY1*YMY1)
            D2=SQRT(XMX2*XMX2+YMY2*YMY2)
            ERR1=ABS(D1+D2-D12(I,L))*D12I(I,L)
            IF (ERR1.LT.ERRMIN) THEN
              IMARK=I
              LMARK=L
              ERRMIN=ERR1
            ENDIF
21          CONTINUE
            IF (I .EQ. NR1STM) THEN
              I = I+1
              GOTO 121
            ENDIF
            HELPP => HELPP%NEXT
          ENDDO
        ENDDO
        IM=IMARK
        LM=LMARK
        IF (ERRMIN.LE.EPS10) GOTO 1000
        IF (LOGX) GOTO 800
        IA=1
        IEP=NP2ND
C
C  NEXT TRY POLOIDAL SURFACES
C  THIS SECTION ALSO: IF LOGY, CHECK ON POLOID. SURFACE IA  (750...)
700     CONTINUE
        DO J=1,4
          HELPP => HELPCUR4(J)%P
          DO WHILE (ASSOCIATED(HELPP))
            I = HELPP%IX
            L = HELPP%IY
221         IF ((I .LT. IAN) .OR. (I .GT. IEN)) GOTO 22
            XMX1=X-XPOL(I,L)
            YMY1=Y-YPOL(I,L)
            XMX4=X-XPOL(I+1,L)
            YMY4=Y-YPOL(I+1,L)
            DX1=SQRT(XMX1*XMX1+YMY1*YMY1)
            DX4=SQRT(XMX4*XMX4+YMY4*YMY4)
            ERR4=ABS(DX1+DX4-D14(I,L))*D14I(I,L)
            IF (ERR4.LT.ERRMIN) THEN
              IMARK=I
              LMARK=L
              ERRMIN=ERR4
            ENDIF
22          CONTINUE
            DO KH=1,NPPLG
              IF (L .EQ. NPOINT(2,KH)-1) THEN
                L = L + 1
                GOTO 221
              ENDIF
            ENDDO
            HELPP => HELPP%NEXT
          ENDDO
        ENDDO
        IM=IMARK
        LM=LMARK
        IF (ERRMIN.LE.EPS10) GOTO 1000
        GOTO 800
C
750     CONTINUE
        DO J=1,4
          HELPP => HELPCUR4(J)%P
          DO WHILE (ASSOCIATED(HELPP))
            I = HELPP%IX
            L = HELPP%IY
            IF (I .GT. NR1STM) GOTO 23
            XMX1=X-XPOL(I,IA)
            YMY1=Y-YPOL(I,IA)
            XMX4=X-XPOL(I+1,IA)
            YMY4=Y-YPOL(I+1,IA)
            DX1=SQRT(XMX1*XMX1+YMY1*YMY1)
            DX4=SQRT(XMX4*XMX4+YMY4*YMY4)
            ERR7=ABS(DX1+DX4-D14(I,IA))*D14I(I,IA)
            IF (ERR7.LT.ERRMIN) THEN
              IMARK=I
              LMARK=IA
              ERRMIN=ERR7
            ENDIF
23          HELPP => HELPP%NEXT
          ENDDO
        ENDDO
        IM=IMARK
        LM=LMARK
        IF (ERRMIN.LE.EPS10) GOTO 1000
C
800     CONTINUE
        IF (INUM.EQ.0.AND.ERRMIN.GT.EPS10) THEN
          CALL MASAGE ('X,Y OUT OF RANGE IN LEARC1                   ')
          CALL MASR2('X,Y             ',X,Y)
          WRITE (6,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
          WRITE (6,*) 'ERRMIN= ',ERRMIN
          WRITE (6,*) 'NPANU,IM,LM= ',NP,IM,LM
          WRITE (6,*) 'IAN,IEN,LOGX,LOGY ',IAN,IEN,LOGX,LOGY
          CALL LEER(1)
        ELSEIF (INUM.GT.1.AND.ERRMIN.GT.EPS10) THEN
          CALL MASAGE ('WARNING FROM LEARC1, INUM.GT.1               ')
          CALL MASR2('X,Y             ',X,Y)
          WRITE (6,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
          WRITE (6,*) 'ERRMIN= ',ERRMIN
          WRITE (6,*) 'NPANU,INUM,IM,LM= ',NP,INUM,IM,LM
          WRITE (6,*) 'IAN,IEN,LOGX,LOGY ',IAN,IEN,LOGX,LOGY
          CALL LEER(1)
        ENDIF
C
1000    CONTINUE
C
        LEARC1=IM
        IPO=LM
C
      ELSEIF (LEVGEO.EQ.2) THEN
C
        LEARC1=0
        IF (LOGX) RETURN
        IF (LOGY) THEN
          IA=1
          IE=NR1STM
        ENDIF
        YQ=Y*Y
        DO 10 J=IA,IE
          I=J+1
          IM=J
          XE=X-EP1(I)
          XEQ=XE*XE
          ATQ=XEQ+YQ/ELLQ(I)
          IF (ATQ.LT.RQ(I)) GOTO 15
10      CONTINUE
        IF (ATQ.LE.RQ(I)+EPS12) GOTO 15
C
        CALL MASAGE ('X,Y OUT OF RANGE IN LEARC1                   ')
        CALL MASR2('X,Y             ',X,Y)
        WRITE (6,*) ATQ,RQ(NR1ST)
        WRITE (6,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
        CALL EXIT_OWN(1)
C
15      CONTINUE
        LEARC1=IM
C
      ELSEIF (LEVGEO.EQ.1) THEN
C
        LEARC1=0
        IF (LOGX) RETURN
        IF (LOGY) THEN
          IA=1
          IE=NR1STM
        ENDIF
        IM=1
        IF (NR1ST.LT.2) GOTO 250
        DO 200 J=IA,IE
          I=J+1
          IM=J
          IF (X.LT.RSURF(I)) GOTO 250
200     CONTINUE
        IF (X.LE.RSURF(I)+EPS12) GOTO 250
C
        CALL MASAGE ('X OUT OF RANGE IN LEARC1                   ')
        CALL MASR2('X,Y             ',X,Y)
        WRITE (6,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
        CALL EXIT_OWN(1)
C
250     CONTINUE
        LEARC1=IM
C
      ELSEIF (LEVGEO.EQ.5) THEN
C
        LEARC1=LEARCT(X,Y,Z)

      ELSEIF (LEVGEO.EQ.6) THEN
C
C  GENERAL GEOMETRY OPTION: PROVIDE CELL NUMBER, GIVEN THE POSITION
C
        LEARC1=LEAUSR(X,Y,Z)
C
      ENDIF
C
      RETURN
      END
