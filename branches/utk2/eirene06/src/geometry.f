C EIRENE06 COMPILATION
C ===== SOURCE: clltst.f
C
      SUBROUTINE CLLTST(*)

      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CGRID
      USE COMPRT

      IMPLICIT NONE

      REAL(DP) :: PHITEST, ZTESTO, PHIT, X0T, Y0T, Z0T, TTT, ZTESTU,
     .          WINK
      INTEGER :: NTEST2, LEARC2, IPOLGT, LEARCT, NTEST0, LEARC1,
     .           LEARCA, LEAUSR, NTEST1

C

      IF (NACELL.GT.0) RETURN
C
C  TEST FOR RADIAL CELL INDICES NRCELL, IPOLG
C
      IF (LEVGEO.LE.4) THEN
C
        NTEST0=LEARC1(X0,Y0,Z0,IPOLGT,1,NR1STM,.FALSE.,.FALSE.,NPANU,
     .                'CLLTST      ')
C
      ELSEIF (LEVGEO.EQ.5) THEN
C
        NTEST0=LEARCT(X0,Y0,Z0)
C
      ELSEIF (LEVGEO.EQ.6) THEN
C
        NTEST0=LEAUSR(X0,Y0,Z0)
C
      ENDIF
C
      IF(NTEST0.NE.NRCELL.AND..NOT.NLSRFX) THEN
        WRITE (iunout,*) 'WRONG CELL-NUMBER IN RADIAL DIRECTION'
        CALL MASJ3 ('NRCELL,NTEST0,NPANU      ',NRCELL,NTEST0,NPANU)
        RETURN 1
      ENDIF
C
C  TEST FOR POLOIDAL CELL INDEX
C
      IF (NLPOL) THEN
        IF (LEVGEO.EQ.1) THEN
          NTEST1=LEARCA(Y0,PSURF,1,NP2ND,1,'CLLTST    ')
        ELSEIF (LEVGEO.EQ.2) THEN
          IF (NLCRC) THEN
            WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
            NTEST1=LEARCA(WINK,PSURF,1,NP2ND,1,'CLLTST    ')
          ELSE
            NTEST1=LEARC2(X0,Y0,NRCELL,NPANU,'CLLTST  ')
          ENDIF
        ELSEIF (LEVGEO.EQ.3) THEN
          NTEST1=IPOLGT
        ELSE
          WRITE (iunout,*) 'ERROR EXIT IN CLLTST, NLPOL ',NPCELL
          CALL EXIT_OWN(1)
        ENDIF
        IF (NTEST1.NE.NPCELL.AND..NOT.NLSRFY) THEN
          WRITE (iunout,*) 'WRONG CELL-NUMBER IN POLOIDAL DIRECTION'
          CALL MASJ3 ('NPCELL,NTEST1,NPANU     ',NPCELL,NTEST1,NPANU)
          RETURN 1
        ENDIF
      ENDIF
C
C  TEST FOR 3RD CELL INDEX
C
      IF (NLTOR) THEN
C  IN CYLINDER APPROXIMATION
        IF (NLTRZ) THEN
          NTEST2=LEARCA(Z0,ZSURF,1,NTTRA,1,'CLLTST      ')
C  IN TOROIDAL APPROXIMATION
        ELSEIF (NLTRA) THEN
          PHIT=PHI
          IF (PHIT.LT.ZSURF(1)) PHIT=PHIT+PI2A
          IF (PHIT.GT.ZSURF(NTTRA)) PHIT=PHIT-PI2A
          PHITEST=ATAN2(Z0,(RMTOR+X0))+ZZONE(NTCELL)
          IF (ABS(PHIT-PHITEST).GT.EPS10) THEN
            WRITE (iunout,*) 'PHI,PHITEST  ',
     .                        PHIT,PHITEST,PHIT-PHITEST,NTCELL
            CALL MASJ1 ('NPANU   ',NPANU)
          ENDIF
          NTEST2=LEARCA(PHIT,ZSURF,1,NTTRA,1,'CLLTST      ')
C
        ENDIF
        IF (NTEST2.NE.NTCELL.AND..NOT.NLSRFZ) THEN
          WRITE (iunout,*) 'WRONG CELL-NUMBER IN TOROIDAL DIRECTION'
          CALL MASJ3 ('NTCELL,NTEST2,NPANU     ',NTCELL,NTEST2,NPANU)
          IF (NLTRA) THEN
            ZTESTO=ZFULL*NTCELL
            ZTESTU=ZFULL*(NTCELL-1)
            CALL MASR3 ('ZTESTU,PHI,ZTESTO=      ',ZTESTU,PHI,ZTESTO)
            RETURN 1
          ENDIF
        ENDIF
C
      ELSEIF (.NOT.NLTOR.AND.NLTRA) THEN
        PHIT=PHI
        IF (PHIT.LT.ZSURF(1)) PHIT=PHIT+PI2A
        IF (PHIT.GT.ZSURF(NTTRA)) PHIT=PHIT-PI2A
        NTEST2=LEARCA(PHIT,ZSURF,1,NTTRA,1,'CLLTST      ')
        PHITEST=ATAN2(Z0,(RMTOR+X0))+ZZONE(NTEST2)
        IF (ABS(PHIT-PHITEST).GT.EPS10) THEN
          WRITE (iunout,*) 'PHI,PHITEST  ',
     .                      PHIT,PHITEST,PHIT-PHITEST,NTEST2
          CALL MASJ1 ('NPANU   ',NPANU)
        ENDIF
C
        X0T=X0+RMTOR
        Z0T=Z0
        TTT=Z0T/(X0T*TANAL)
        IF (ABS(TTT).GT.1.+EPS10) THEN
          WRITE (iunout,*) 'WRONG CO-ORDINATES IN TOROIDAL DIRECTION'
          CALL MASJ1 ('NPANU   ',NPANU)
          CALL MASR3 ('X0,Z0,TTT                ',X0,Z0,TTT)
          RETURN 1
        ENDIF
C
        IF (NTEST2.NE.IPERID.AND..NOT.NLSRFZ) THEN
          WRITE (iunout,*) 'WRONG CELL-NUMBER IN TOROIDAL DIRECTION'
          CALL MASJ3 ('IPERID,NTEST2,NPANU     ',IPERID,NTEST2,NPANU)
        ENDIF
      ENDIF
C
      RETURN
      END
C ===== SOURCE: fzrtor.f
C
C
      SUBROUTINE FZRTOR(X,Z,NOLD,XNEW,PH,NNEW,LTEST,NTEST)
C
C  INPUT :
C    X     : X CO-ORDINATE
C    Z     : Z CO-ORDINATE
C    NOLD  : X AND Z ARE GIVEN IN LOCAL SYSTEM NOLD
C  OUTPUT:
C    PH   : TOROIDAL ANGLE (0<=PH<2PI)
C    NNEW  : TOROIDAL ZONE NUMBER (1 <= NNEW <= NTTRAM)
C    XNEW  : X CO-ORDINATE IN LOCAL SYSTEM OF CELL NNEW
C  IF LTEST: TEST IF NNEW=NTEST, NTEST IS INPUT (EXPECTED VALUE OF NNEW)
C
C  FROM LOCAL CO-ORDINATES X,Z IN THE TOROIDAL CELL NUMBER NOLD
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGRID
      USE COMPRT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X, Z
      REAL(DP), INTENT(OUT) :: XNEW, PH
      INTEGER, INTENT(IN) :: NOLD, NTEST
      INTEGER, INTENT(OUT) :: NNEW
      LOGICAL, INTENT(IN) :: LTEST
      REAL(DP) :: Z1, ZOLD, XOLD, WLOC, RR
      INTEGER :: LEARCA
C
      XOLD=X+RMTOR
      ZOLD=Z
C  TOROIDAL ANGLE PH AND NEW TOROIDAL CELL NUMBER
      Z1=ZSURF(1)
      PH=MOD(ZZONE(NOLD)+ATAN2(ZOLD,XOLD)+PI2A-Z1,PI2A)+Z1
      NNEW=LEARCA(PH,ZSURF,1,NTTRA,1,'FZRTOR ')
      NNEW=MIN(NTTRAM,NNEW)
C  RADIAL CO-ORDINATE
      RR=SQRT(XOLD*XOLD+ZOLD*ZOLD)
C  LOCAL CO-ORDINATES IN CELL NNEW
      WLOC=PH-ZZONE(NNEW)
      XNEW=RR*COS(WLOC)
      IF (LTEST) THEN
        IF (NTEST.NE.NNEW) THEN
C  POSITION ON (NEAR) TOROIDAL SURFACE?
          IF (ABS(PH-ZSURF(NTEST)).GT.EPS12.AND.
     .        ABS(PH-ZSURF(NTEST+1)).GT.EPS12) THEN
            WRITE (iunout,*) 
     .       'ERROR IN FZRTOR: WRONG TOROIDAL CELL INDEX'
            WRITE (iunout,*) 'XOLD,ZOLD,XNEW,PH  ',XOLD,ZOLD,XNEW,PH
            WRITE (iunout,*) 'NNEW,NTST ',NNEW,NTEST
            WRITE (iunout,*) 'ABS(PH-ZSURF(NTST)),ABS(PH-ZSURF(NTST+1))'
            WRITE (iunout,*) 
     .             ABS(PH-ZSURF(NTEST)),ABS(PH-ZSURF(NTEST+1))
            WEIGHT=0.
            LGPART=.FALSE.
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
C ===== SOURCE: fzrtra.f
!pb  05.10.06 : intent of arguments X and Z corrected
C
C
      SUBROUTINE FZRTRA(X,Z,PH,NNEW)

      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGRID
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: X, Z
      REAL(DP), INTENT(INOUT) :: PH
      INTEGER, INTENT(OUT) :: NNEW
      REAL(DP) :: Z1, X01,XX
      INTEGER :: LEARCA
C
C
C  FIND X,Z, NNEW,   FROM X,PH   (X=X??)
      Z1=ZSURF(1)
      PH=MOD(PH+PI2A-Z1,PI2A)+Z1
      NNEW=LEARCA(PH,ZSURF,1,NTTRA,1,'FZRTRA ')
      IF (NNEW.LE.0.OR.NNEW.GT.NTTRAM) THEN
        WRITE (iunout,*) 'NT OUT OF RANGE IN FZRTRA '
        WRITE (iunout,*) PH,ZHALF,NNEW
        CALL EXIT_OWN(1)
      ENDIF
      X01=X+RMTOR
      CALL FZRTRI(X,Z,NNEW,X01,PH,NNEW)
      RETURN
      END
C ===== SOURCE: fzrtri.f
C
C
      SUBROUTINE FZRTRI(X0,Z0,NNEW,XOLD,PH,NOLD)
C
C  CALCULATE THE LOCAL CO-ORDINATES X0,Z0 IN ZONE NNEW
C
C  FROM THE CELL NUMBER NOLD,
C       THE X CO-ORDINATE XOLD (=X+RMTOR) IN THAT CELL
C       THE TOROIDAL ANGLE PH
C
C  INPUT : XOLD,PH,NOLD, NNEW
C  OUTPUT: X0,Z0
C
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGRID

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XOLD, PH
      REAL(DP), INTENT(OUT) :: X0, Z0
      INTEGER, INTENT(IN) :: NOLD, NNEW
      REAL(DP) :: WLOC1, WLOC2, WTH, Z1, RR, XNEW

      Z1=ZSURF(1)
      WTH=MOD(PH+PI2A-Z1,PI2A)+Z1
      WLOC1=WTH-ZZONE(NOLD)
      RR=XOLD/COS(WLOC1)
C
      WLOC2=WTH-ZZONE(NNEW)
      XNEW=RR*COS(WLOC2)
      X0=XNEW-RMTOR
      Z0=TAN(WLOC2)*XNEW
      RETURN
      END
C ===== SOURCE: learc1.f
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
      USE COMPRT, ONLY: IUNOUT
      USE CGRID
      USE CGEOM
      USE CCONA
      USE CPOLYG
      USE CLOGAU
      USE CTRIG

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X, Y, Z
      INTEGER, INTENT(IN) :: IAN, IEN, NP
      INTEGER, INTENT(INOUT) :: IPO
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
                WRITE (iunout,*) 'OBSCURE CELL DETECTED, ITRII= ',I
                WRITE (iunout,*) 'OPPOSITE ORIENTATION IN TRIANGLE '
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
        IM = 0
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
          WRITE (iunout,*) 'NO TRIANGLE FOUND IN LEARC1 FOR '
          WRITE (iunout,*) 'X = ',X,' Y = ',Y
          WRITE (iunout,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
          WRITE (iunout,*) 'NPANU,IM= ',NP,IM
        ELSEIF (INUM.GT.1) THEN
          CALL MASAGE ('WARNING FROM LEARC1, INUM.GT.1               ')
          CALL MASR2('X,Y             ',X,Y)
          WRITE (iunout,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
          WRITE (iunout,*) 'NPANU,INUM,IM= ',NP,INUM,IM
          WRITE (iunout,*) 'IAN,IEN,LOGX,LOGY ',IAN,IEN,LOGX,LOGY
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
                WRITE (iunout,*) 'OBSCURE CELL DETECTED, IR,IP= ',I,L
                WRITE (iunout,*) 
     .            'PLASMA FIELD LIKELY TO BE CORRUPTED HERE '
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
              WRITE (iunout,*) 'OPPOSITE GRID ORIENTATION '
              WRITE (iunout,*) 'FUNCTION LEARC1 NEEDS TO BE EXTENDED'
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
          WRITE (iunout,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
          WRITE (iunout,*) 'ERRMIN= ',ERRMIN
          WRITE (iunout,*) 'NPANU,IM,LM= ',NP,IM,LM
          WRITE (iunout,*) 'IAN,IEN,LOGX,LOGY ',IAN,IEN,LOGX,LOGY
          CALL LEER(1)
        ELSEIF (INUM.GT.1.AND.ERRMIN.GT.EPS10) THEN
          CALL MASAGE ('WARNING FROM LEARC1, INUM.GT.1               ')
          CALL MASR2('X,Y             ',X,Y)
          WRITE (iunout,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
          WRITE (iunout,*) 'ERRMIN= ',ERRMIN
          WRITE (iunout,*) 'NPANU,INUM,IM,LM= ',NP,INUM,IM,LM
          WRITE (iunout,*) 'IAN,IEN,LOGX,LOGY ',IAN,IEN,LOGX,LOGY
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
        WRITE (iunout,*) ATQ,RQ(NR1ST)
        WRITE (iunout,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
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
        WRITE (iunout,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
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
C ===== SOURCE: learc2.f
C
C
      FUNCTION LEARC2(X,Y,NR,NP,TEXT)
C
C  THIS SUBROUTINE FINDS THE POLYGON INDEX "IPOLG"
C  ASSUMING THAT X,Y IS IN THE RADIAL ZONE NR
C
      USE PRECISION
      USE PARMMOD
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CPOLYG
      USE CGRID
      USE CGEOM

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X, Y
      INTEGER, INTENT(IN) :: NR, NP
      CHARACTER(LEN=*), INTENT(IN) :: TEXT

      REAL(DP) :: ERR1(N2NDPLG), ERR2(N2NDPLG), ERR3(N2NDPLG),
     .          ERR4(N2NDPLG), ERR5(N2NDPLG), ERR6(N2NDPLG)
      REAL(DP) :: HELPN, DWYN, WY1N, XMX3, ERRMIN, UXN, UYN, TXN, TYN,
     .          WX1N, VY1N, DETN, YMY3, X2N, X4N,
     .          DX1, DX2, DX3, DX4, XMX4, YMY4, VX1, VX2, WY2,
     .          WX2, VY2, Y3N, X3N, Y1N, Y2N, Y4N, X1N, DET1,
     .          YMY2, XMX2, DET2, YMY1, XMX1
      INTEGER :: IM, LM, LMARK, IMARK, L, K, LEARC2, N, IFIRST, INUM
      REAL(DP), ALLOCATABLE, SAVE ::
     .          X1(:,:),Y1(:,:),X2(:,:),Y2(:,:),X3(:,:),Y3(:,:),
     .          X4(:,:),Y4(:,:),TX(:,:),TY(:,:),
     .          UX(:,:),UY(:,:),DET(:,:),
     .          VY1(:,:),WY1(:,:),WX1(:,:),DWY(:,:),HELP(:,:),
     .          D12(:,:),D14(:,:),D32(:,:),D34(:,:)
!pb      SAVE
      DATA IFIRST /0/
C
      IF (IFIRST .EQ. 0) THEN
        IFIRST = 1
        ALLOCATE (X1(N1STS,N2NDS))
        ALLOCATE (Y1(N1STS,N2NDS))
        ALLOCATE (X2(N1STS,N2NDS))
        ALLOCATE (Y2(N1STS,N2NDS))
        ALLOCATE (X3(N1STS,N2NDS))
        ALLOCATE (Y3(N1STS,N2NDS))
        ALLOCATE (X4(N1STS,N2NDS))
        ALLOCATE (Y4(N1STS,N2NDS))
        ALLOCATE (TX(N1STS,N2NDS))
        ALLOCATE (TY(N1STS,N2NDS))
        ALLOCATE (UX(N1STS,N2NDS))
        ALLOCATE (UY(N1STS,N2NDS))
        ALLOCATE (DET(N1STS,N2NDS))
        ALLOCATE (VY1(N1STS,N2NDS))
        ALLOCATE (WY1(N1STS,N2NDS))
        ALLOCATE (WX1(N1STS,N2NDS))
        ALLOCATE (DWY(N1STS,N2NDS))
        ALLOCATE (HELP(N1STS,N2NDS))
        ALLOCATE (D12(N1STS,N2NDS))
        ALLOCATE (D14(N1STS,N2NDS))
        ALLOCATE (D32(N1STS,N2NDS))
        ALLOCATE (D34(N1STS,N2NDS))
        DO 10 N=1,NR1STM
          DO 10 K=1,NPPLG
            DO 10 L=NPOINT(1,K),NPOINT(2,K)-1
              X1(N,L)=XPOL(N,L)
              Y1(N,L)=YPOL(N,L)
              X2(N,L)=XPOL(N,L+1)
              Y2(N,L)=YPOL(N,L+1)
              X3(N,L)=XPOL(N+1,L+1)
              Y3(N,L)=YPOL(N+1,L+1)
              X4(N,L)=XPOL(N+1,L)
              Y4(N,L)=YPOL(N+1,L)
              TX(N,L)=X1(N,L)-X2(N,L)
              TY(N,L)=Y1(N,L)-Y2(N,L)
              UX(N,L)=X3(N,L)-X2(N,L)
              UY(N,L)=Y3(N,L)-Y2(N,L)
              DET(N,L)=1./(TX(N,L)*UY(N,L)-TY(N,L)*UX(N,L)-EPS60)
C
              VX1=X4(N,L)-X1(N,L)
              VY1(N,L)=Y4(N,L)-Y1(N,L)
              WX1(N,L)=X3(N,L)-X1(N,L)
              WY1(N,L)=Y3(N,L)-Y1(N,L)
              DWY(N,L)=1./(WY1(N,L)+EPS60)
              HELP(N,L)=1./(VX1*WY1(N,L)-VY1(N,L)*WX1(N,L)-EPS60)
C
              VX2=X3(N,L)-X4(N,L)
              VY2=Y3(N,L)-Y4(N,L)
              WX2=X1(N,L)-X4(N,L)
              WY2=Y1(N,L)-Y4(N,L)
              D12(N,L)=SQRT(TX(N,L)*TX(N,L)+TY(N,L)*TY(N,L))
              D32(N,L)=SQRT(UX(N,L)*UX(N,L)+UY(N,L)*UY(N,L))
              D34(N,L)=SQRT(VX2*VX2+VY2*VY2)
              D14(N,L)=SQRT(WX2*WX2+WY2*WY2)
10      CONTINUE
      ENDIF
C
C  END OF IFIRST LOOP
C
      INUM=0
C
      DO 100 K=1,NPPLG
      DO 101 L=NPOINT(1,K),NPOINT(2,K)-1
        XMX2=X-X2(NR,L)
        YMY2=Y-Y2(NR,L)
        DET1=XMX2*UY(NR,L)-YMY2*UX(NR,L)
        DET2=TX(NR,L)*YMY2-TY(NR,L)*XMX2
        ERR1(L)=DET1*DET(NR,L)
        ERR2(L)=DET2*DET(NR,L)
C
        XMX1=X-X1(NR,L)
        YMY1=Y-Y1(NR,L)
        ERR3(L)=(XMX1*WY1(NR,L)-YMY1*WX1(NR,L))*HELP(NR,L)
        ERR4(L)=(YMY1-VY1(NR,L)*ERR3(L))*DWY(NR,L)
        ERR5(L)=ERR1(L)+ERR2(L)
        ERR6(L)=ERR3(L)+ERR4(L)
101   CONTINUE
      IF (LEVGEO .EQ. 2) THEN
        DO 102 L=NPOINT(1,K),NPOINT(2,K)-1
102       ERR6(L)=0.
      ENDIF
      DO 100 L=NPOINT(1,K),NPOINT(2,K)-1
        IF (ERR4(L).GT.1.D30) THEN
          X2N=XPOL(NR,L)
          Y2N=YPOL(NR,L)
          X3N=XPOL(NR,L+1)
          Y3N=YPOL(NR,L+1)
          X4N=XPOL(NR+1,L+1)
          Y4N=YPOL(NR+1,L+1)
          X1N=XPOL(NR+1,L)
          Y1N=YPOL(NR+1,L)
          TXN=X1N-X2N
          TYN=Y1N-Y2N
          UXN=X3N-X2N
          UYN=Y3N-Y2N
          DETN=1./(TXN*UYN-TYN*UXN-EPS60)

          VX1=X4N-X1N
          VY1N=Y4N-Y1N
          WX1N=X3N-X1N
          WY1N=Y3N-Y1N
          DWYN=1./(WY1N+EPS60)
          HELPN=1./(VX1*WY1N-VY1N*WX1N-EPS60)

          XMX2=X-X2N
          YMY2=Y-Y2N
          DET1=XMX2*UYN-YMY2*UXN
          DET2=TXN*YMY2-TYN*XMX2
          ERR1(L)=DET1*DETN
          ERR2(L)=DET2*DETN
C
          XMX1=X-X1N
          YMY1=Y-Y1N
          ERR3(L)=(XMX1*WY1N-YMY1*WX1N)*HELPN
          ERR4(L)=(YMY1-VY1N*ERR3(L))*DWYN
          ERR5(L)=ERR1(L)+ERR2(L)
          ERR6(L)=ERR3(L)+ERR4(L)
        ENDIF
        IF (ERR5(L).LT.1..AND.ERR1(L).GT.0..AND.ERR2(L).GT.0.) THEN
          INUM=INUM+1
          IM=NR
          LM=L
        ELSEIF (ERR6(L).LT.1..AND.ERR3(L).GE.0..AND.ERR4(L).GE.0.) THEN
          INUM=INUM+1
          IM=NR
          LM=L
        ENDIF
100   CONTINUE
C
      ERRMIN=1.D30
      IF (INUM.NE.1) THEN
C  CHECK FOR NEAREST BOUNDARY, BECAUSE NO VALID CELL INDEX FOUND
        DO 110 K=1,NPPLG
        DO 111 L=NPOINT(1,K),NPOINT(2,K)-1
          XMX1=X-X1(NR,L)
          YMY1=Y-Y1(NR,L)
          XMX2=X-X2(NR,L)
          YMY2=Y-Y2(NR,L)
          XMX3=X-X3(NR,L)
          YMY3=Y-Y3(NR,L)
          XMX4=X-X4(NR,L)
          YMY4=Y-Y4(NR,L)
          DX1=SQRT(XMX1*XMX1+YMY1*YMY1)
          DX2=SQRT(XMX2*XMX2+YMY2*YMY2)
          DX3=SQRT(XMX3*XMX3+YMY3*YMY3)
          DX4=SQRT(XMX4*XMX4+YMY4*YMY4)
          ERR1(L)=ABS(DX1+DX2-D12(NR,L))
          ERR2(L)=ABS(DX2+DX3-D32(NR,L))
          ERR3(L)=ABS(DX3+DX4-D34(NR,L))
          ERR4(L)=ABS(DX1+DX4-D14(NR,L))
111     CONTINUE
        IF (LEVGEO .EQ. 2) THEN
!PB          ERR1(L)=ERRMIN
!PB          ERR3(L)=ERRMIN
          ERR1=ERRMIN
          ERR3=ERRMIN
        ENDIF
        DO 110 L=NPOINT(1,K),NPOINT(2,K)-1
          IF (ERR1(L).LT.ERRMIN) THEN
            IMARK=NR
            LMARK=L
            ERRMIN=ERR1(L)
          ENDIF
          IF (ERR2(L).LT.ERRMIN) THEN
            IMARK=NR
            LMARK=L+1
            ERRMIN=ERR2(L)
          ENDIF
          IF (ERR3(L).LT.ERRMIN) THEN
            IMARK=NR+1
            LMARK=L
            ERRMIN=ERR3(L)
          ENDIF
          IF (ERR4(L).LT.ERRMIN) THEN
            IMARK=NR
            LMARK=L
            ERRMIN=ERR4(L)
          ENDIF
110     CONTINUE
        IM=IMARK
        LM=LMARK
      ENDIF
C
      IF (INUM.EQ.0.AND.ERRMIN.GT.EPS5) THEN
        CALL MASAGE ('X,Y OUT OF RANGE IN LEARC2                   ')
        CALL MASR2('X,Y             ',X,Y)
        WRITE (iunout,*) 'LEARC2 CALLED FROM SUBR. ',TEXT
        WRITE (iunout,*) 'ERRMIN= ',ERRMIN
        WRITE (iunout,*) 'NPANU,IM,LM= ',NP,IM,LM
      ELSEIF (INUM.GT.1.AND.ERRMIN.GT.EPS5) THEN
        CALL MASAGE ('WARNING FROM LEARC2, INUM.GT.1               ')
        CALL MASR2('X,Y             ',X,Y)
        WRITE (iunout,*) 'LEARC2 CALLED FROM SUBR. ',TEXT
        WRITE (iunout,*) 'ERRMIN= ',ERRMIN
        WRITE (iunout,*) 'NPANU,INUM,IM,LM= ',NP,INUM,IM,LM
      ENDIF
C
      LEARC2=LM
C
      RETURN
      END
C ===== SOURCE: rotadd.f
C
C
      SUBROUTINE ROTADD(A,AI,ILINI,ILEND)
C
C  TRANSFORM CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES
C  CO-ORDINATES OF A POINT IN OLD SYSTEM (X,Y,Z)
C  CO-ORDINATES OF A POINT IN NEW SYSTEM (X',Y',Z')
C
C    X'         X     X          X'                           -1  T
C    Y' =   A * Y  ;  Y  =  AI * Y' ; SUBROUTINE ASSUMES: AI=A  =A
C    Z'         Z     Z          Z'
C
      USE PRECISION
      USE PARMMOD
      USE COMPRT, ONLY: IUNOUT
      USE CADGEO
      USE CTRCEI

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A(3,3),AI(3,3)
      INTEGER,  INTENT(IN) :: ILINI, ILEND
      REAL(DP) :: P1S, P2S, P3S, A3S, A4S, A5S, A6S, A7S, A8S, A9S, S,
     .            DEL, DETER, T, A1S, A2S
      INTEGER  :: I, J
C
C
      DO 100 I=ILINI,ILEND
C
C   CHANGE COEFFICIENTS FOR ALGEBRAIC EQUATION, TRY TO KEEP INVARIANTS
          S=A4LM(I)+A5LM(I)+A6LM(I)
          DEL=DETER(A4LM(I),A7LM(I)/2._DP,A9LM(I)/2._DP,A7LM(I)/2._DP,
     .              A5LM(I),A8LM(I)/2._DP,A9LM(I)/2._DP,A8LM(I)/2._DP,
     .              A6LM(I))
          T=A5LM(I)*A6LM(I)+A6LM(I)*A4LM(I)+A4LM(I)*A5LM(I)-
     .     (A9LM(I)**2+A8LM(I)**2+A7LM(I)**2)*0.25
C         IF (TRCPLT) THEN
C           WRITE (iunout,*) 'INVARIANTS FROM ROTADD: I= ',I
C           WRITE (iunout,*) 'BEFORE: S,DEL,T= ',S,DEL,T
C           WRITE (iunout,*) ' A0...A9 ',
C    .                  A0LM(I),A1LM(I),A2LM(I),A3LM(I),
C    .                  A4LM(I),A5LM(I),A6LM(I),A7LM(I),A8LM(I),A9LM(I)
C         ENDIF
C
C         A0LM(I)=A0LM(I)
          A1S=A1LM(I)*AI(1,1)+A2LM(I)*AI(2,1)+A3LM(I)*AI(3,1)
          A2S=A1LM(I)*AI(1,2)+A2LM(I)*AI(2,2)+A3LM(I)*AI(3,2)
          A3S=A1LM(I)*AI(1,3)+A2LM(I)*AI(2,3)+A3LM(I)*AI(3,3)
          A1LM(I)=A1S
          A2LM(I)=A2S
          A3LM(I)=A3S
C
          A4S=(A4LM(I)*AI(1,1)+A7LM(I)*AI(2,1))*AI(1,1)+
     .        (A5LM(I)*AI(2,1)+A9LM(I)*AI(3,1))*AI(2,1)+
     .        (A6LM(I)*AI(3,1)+A8LM(I)*AI(1,1))*AI(3,1)
          A5S=(A4LM(I)*AI(1,2)+A7LM(I)*AI(2,2))*AI(1,2)+
     .        (A5LM(I)*AI(2,2)+A9LM(I)*AI(3,2))*AI(2,2)+
     .        (A6LM(I)*AI(3,2)+A8LM(I)*AI(1,2))*AI(3,2)
          A6S=(A4LM(I)*AI(1,3)+A7LM(I)*AI(2,3))*AI(1,3)+
     .        (A5LM(I)*AI(2,3)+A9LM(I)*AI(3,3))*AI(2,3)+
     .        (A6LM(I)*AI(3,3)+A8LM(I)*AI(1,3))*AI(3,3)
C         A6S=S-A4S-A5S
C
          A7S=A4LM(I)*AI(1,1)*AI(1,2)*2+
     .        A5LM(I)*AI(2,1)*AI(2,2)*2+
     .        A6LM(I)*AI(3,1)*AI(3,2)*2+
     .        A7LM(I)*(AI(1,1)*AI(2,2)+AI(1,2)*AI(2,1))+
     .        A8LM(I)*(AI(1,1)*AI(3,2)+AI(1,2)*AI(3,1))+
     .        A9LM(I)*(AI(2,1)*AI(3,2)+AI(2,2)*AI(3,1))
          A8S=A4LM(I)*AI(1,1)*AI(1,3)*2+
     .        A5LM(I)*AI(2,1)*AI(2,3)*2+
     .        A6LM(I)*AI(3,1)*AI(3,3)*2+
     .        A7LM(I)*(AI(1,1)*AI(2,3)+AI(1,3)*AI(2,1))+
     .        A8LM(I)*(AI(1,1)*AI(3,3)+AI(1,3)*AI(3,1))+
     .        A9LM(I)*(AI(2,1)*AI(3,3)+AI(2,3)*AI(3,1))
          A9S=A4LM(I)*AI(1,2)*AI(1,3)*2+
     .        A5LM(I)*AI(2,2)*AI(2,3)*2+
     .        A6LM(I)*AI(3,2)*AI(3,3)*2+
     .        A7LM(I)*(AI(2,1)*AI(2,3)+AI(1,3)*AI(2,2))+
     .        A8LM(I)*(AI(1,2)*AI(3,3)+AI(1,3)*AI(3,2))+
     .        A9LM(I)*(AI(2,2)*AI(3,3)+AI(2,3)*AI(3,2))
          A4LM(I)=A4S
          A5LM(I)=A5S
          A6LM(I)=A6S
          A7LM(I)=A7S
          A8LM(I)=A8S
          A9LM(I)=A9S
C
          S=A4LM(I)+A5LM(I)+A6LM(I)
          DEL=DETER(A4LM(I),A7LM(I)/2._DP,A9LM(I)/2._DP,A7LM(I)/2._DP,
     .              A5LM(I),A8LM(I)/2._DP,A9LM(I)/2._DP,A8LM(I)/2._DP,
     .              A6LM(I))
          T=A5LM(I)*A6LM(I)+A6LM(I)*A4LM(I)+A4LM(I)*A5LM(I)-
     .     (A9LM(I)**2+A8LM(I)**2+A7LM(I)**2)*0.25
C         WRITE (iunout,*) 'AFTER:  S,DEL,T= ',S,DEL,T
C
1         IF (RLB(I).LT.0.) THEN
C   CHANGE COEFFICIENTS IN LINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 10 J=1,ILIN(I)
C             ALIMS(J,I)=ALIMS(J,I)
              A1S=XLIMS(J,I)*AI(1,1)+YLIMS(J,I)*AI(2,1)+
     .            ZLIMS(J,I)*AI(3,1)
              A2S=XLIMS(J,I)*AI(1,2)+YLIMS(J,I)*AI(2,2)+
     .            ZLIMS(J,I)*AI(3,2)
              A3S=XLIMS(J,I)*AI(1,3)+YLIMS(J,I)*AI(2,3)+
     .            ZLIMS(J,I)*AI(3,3)
              XLIMS(J,I)=A1S
              YLIMS(J,I)=A2S
              ZLIMS(J,I)=A3S
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
C             ALIMS0(J,I)=ALIMS0(J,I)
              A1S=XLIMS1(J,I)*AI(1,1)+YLIMS1(J,I)*AI(2,1)+
     .            ZLIMS1(J,I)*AI(3,1)
              A2S=XLIMS1(J,I)*AI(1,2)+YLIMS1(J,I)*AI(2,2)+
     .            ZLIMS1(J,I)*AI(3,2)
              A3S=XLIMS1(J,I)*AI(1,3)+YLIMS1(J,I)*AI(2,3)+
     .            ZLIMS1(J,I)*AI(3,3)
              XLIMS1(J,I)=A1S
              YLIMS1(J,I)=A2S
              ZLIMS1(J,I)=A3S
C
              S=XLIMS2(J,I)+YLIMS2(J,I)+ZLIMS2(J,I)
              A4S=XLIMS2(J,I)*AI(1,1)**2+YLIMS2(J,I)*AI(2,1)**2+
     .            ZLIMS2(J,I)*AI(3,1)**2+
     .            XLIMS3(J,I)*AI(1,1)*AI(2,1)+
     .            YLIMS3(J,I)*AI(1,1)*AI(3,1)+
     .            ZLIMS3(J,I)*AI(2,1)*AI(3,1)
              A5S=XLIMS2(J,I)*AI(1,2)**2+YLIMS2(J,I)*AI(2,2)**2+
     .            ZLIMS2(J,I)*AI(3,2)**2+
     .            XLIMS3(J,I)*AI(1,2)*AI(2,2)+
     .            YLIMS3(J,I)*AI(1,2)*AI(3,2)+
     .            ZLIMS3(J,I)*AI(2,2)*AI(3,2)
              A6S=S-A4S-A5S
C
              A7S=XLIMS2(J,I)*(AI(1,1)*AI(1,2)+AI(1,1)*AI(1,2))+
     .            YLIMS2(J,I)*(AI(2,1)*AI(2,2)+AI(2,1)*AI(2,2))+
     .            ZLIMS2(J,I)*(AI(3,1)*AI(3,2)+AI(3,1)*AI(3,2))+
     .            XLIMS3(J,I)*(AI(1,1)*AI(2,2)+AI(1,2)*AI(2,1))+
     .            YLIMS3(J,I)*(AI(1,1)*AI(3,2)+AI(1,2)*AI(3,1))+
     .            ZLIMS3(J,I)*(AI(2,1)*AI(3,2)+AI(2,2)*AI(3,1))
              A8S=XLIMS2(J,I)*(AI(1,1)*AI(1,3)+AI(1,1)*AI(1,3))+
     .            YLIMS2(J,I)*(AI(2,1)*AI(2,3)+AI(2,1)*AI(2,3))+
     .            ZLIMS2(J,I)*(AI(3,1)*AI(3,3)+AI(3,1)*AI(3,3))+
     .            XLIMS3(J,I)*(AI(1,1)*AI(2,3)+AI(1,3)*AI(2,1))+
     .            YLIMS3(J,I)*(AI(1,1)*AI(3,3)+AI(1,3)*AI(3,1))+
     .            ZLIMS3(J,I)*(AI(2,1)*AI(3,3)+AI(2,3)*AI(3,1))
              A9S=XLIMS2(J,I)*(AI(1,2)*AI(1,3)+AI(1,2)*AI(1,3))+
     .            YLIMS2(J,I)*(AI(2,2)*AI(2,3)+AI(2,2)*AI(2,3))+
     .            ZLIMS2(J,I)*(AI(3,2)*AI(3,3)+AI(3,2)*AI(3,3))+
     .            XLIMS3(J,I)*(AI(2,1)*AI(2,3)+AI(1,3)*AI(2,2))+
     .            YLIMS3(J,I)*(AI(1,2)*AI(3,3)+AI(1,3)*AI(3,2))+
     .            ZLIMS3(J,I)*(AI(2,2)*AI(3,3)+AI(2,3)*AI(3,2))
              XLIMS2(J,I)=A4S
              YLIMS2(J,I)=A5S
              ZLIMS2(J,I)=A6S
              XLIMS3(J,I)=A7S
              YLIMS3(J,I)=A8S
              ZLIMS3(J,I)=A9S
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER: CHANGE TO RLB(I)=-6 OPTION AND DEFINE 6 LINEAR
C                     INEQUALITIES
            ALIMS(1,I)=XLIMS1(1,I)
            XLIMS(1,I)=-1.
            YLIMS(1,I)=0.
            ZLIMS(1,I)=0.
            ALIMS(2,I)=-XLIMS2(1,I)
            XLIMS(2,I)=1.
            YLIMS(2,I)=0.
            ZLIMS(2,I)=0.
            ALIMS(3,I)=YLIMS1(1,I)
            XLIMS(3,I)=0.
            YLIMS(3,I)=-1.
            ZLIMS(3,I)=0.
            ALIMS(4,I)=-YLIMS2(1,I)
            XLIMS(4,I)=0.
            YLIMS(4,I)=1.
            ZLIMS(4,I)=0.
            ALIMS(5,I)=ZLIMS1(1,I)
            XLIMS(5,I)=0.
            YLIMS(5,I)=0.
            ZLIMS(5,I)=-1.
            ALIMS(6,I)=-ZLIMS2(1,I)
            XLIMS(6,I)=0.
            YLIMS(6,I)=0.
            ZLIMS(6,I)=1.
            IF (RLB(I).EQ.1.5) THEN
              WRITE (iunout,*) 'ROTADD: TO BE WRITTEN, RLB=1.5'
            ENDIF
            RLB(I)=-6.
            ILIN(I)=6
            ISCN(I)=0
            GOTO 1
        ELSE
C
C   POINT OPTIONS
          P1S=P1(1,I)*A(1,1)+P1(2,I)*A(1,2)+P1(3,I)*A(1,3)
          P2S=P1(1,I)*A(2,1)+P1(2,I)*A(2,2)+P1(3,I)*A(2,3)
          P3S=P1(1,I)*A(3,1)+P1(2,I)*A(3,2)+P1(3,I)*A(3,3)
          P1(1,I)=P1S
          P1(2,I)=P2S
          P1(3,I)=P3S
          P1S=P2(1,I)*A(1,1)+P2(2,I)*A(1,2)+P2(3,I)*A(1,3)
          P2S=P2(1,I)*A(2,1)+P2(2,I)*A(2,2)+P2(3,I)*A(2,3)
          P3S=P2(1,I)*A(3,1)+P2(2,I)*A(3,2)+P2(3,I)*A(3,3)
          P2(1,I)=P1S
          P2(2,I)=P2S
          P2(3,I)=P3S
          IF (RLB(I).GE.3.) THEN
            P1S=P3(1,I)*A(1,1)+P3(2,I)*A(1,2)+P3(3,I)*A(1,3)
            P2S=P3(1,I)*A(2,1)+P3(2,I)*A(2,2)+P3(3,I)*A(2,3)
            P3S=P3(1,I)*A(3,1)+P3(2,I)*A(3,2)+P3(3,I)*A(3,3)
            P3(1,I)=P1S
            P3(2,I)=P2S
            P3(3,I)=P3S
          ENDIF
          IF (RLB(I).GE.4.) THEN
            P1S=P4(1,I)*A(1,1)+P4(2,I)*A(1,2)+P4(3,I)*A(1,3)
            P2S=P4(1,I)*A(2,1)+P4(2,I)*A(2,2)+P4(3,I)*A(2,3)
            P3S=P4(1,I)*A(3,1)+P4(2,I)*A(3,2)+P4(3,I)*A(3,3)
            P4(1,I)=P1S
            P4(2,I)=P2S
            P4(3,I)=P3S
          ENDIF
          IF (RLB(I).GE.5.) THEN
            P1S=P5(1,I)*A(1,1)+P5(2,I)*A(1,2)+P5(3,I)*A(1,3)
            P2S=P5(1,I)*A(2,1)+P5(2,I)*A(2,2)+P5(3,I)*A(2,3)
            P3S=P5(1,I)*A(3,1)+P5(2,I)*A(3,2)+P5(3,I)*A(3,3)
            P5(1,I)=P1S
            P5(2,I)=P2S
            P5(3,I)=P3S
          ENDIF
          IF (RLB(I).GE.6.) THEN
            P1S=P6(1,I)*A(1,1)+P6(2,I)*A(1,2)+P6(3,I)*A(1,3)
            P2S=P6(1,I)*A(2,1)+P6(2,I)*A(2,2)+P6(3,I)*A(2,3)
            P3S=P6(1,I)*A(3,1)+P6(2,I)*A(3,2)+P6(3,I)*A(3,3)
            P6(1,I)=P1S
            P6(2,I)=P2S
            P6(3,I)=P3S
          ENDIF
        ENDIF
C
100   CONTINUE
C
      RETURN
      END
C ===== SOURCE: sneigh.f
C
C
      SUBROUTINE SNEIGH
C
C  DETERMINE THE INDICES OF THE FOUR NEIGHBORING CELLS OF
C  AN EIRENE CELL    : NGHPOL
C  AN EIRENE SURFCASE: NGHPLS
C  SIDE NUMBERING:
C
C       (IR+1,IP+1)          (IR,IP+1)
C                       2
C                 +-----------+
C                 |           |
C                 |           |
C               3 |           | 1
C                 |           |
C                 |           |
C                 +-----------+
C                       4
C         (IR+1,IP)          (IR,IP)
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE module_avltree

      IMPLICIT NONE

      INTEGER :: IR, IP, IPART, JP, K, IC, IN
      TYPE(CELL_ELEM), POINTER :: CUR
      type(TAVLTree), pointer :: baum
      logical :: inserted

      DO IR=1,NR1ST
        DO IP=1,NP2ND
          NGHPOL(1,IR,IP)=0
          NGHPOL(2,IR,IP)=0
          NGHPOL(3,IR,IP)=0
          NGHPOL(4,IR,IP)=0
          NGHPLS(1,IR,IP)=0
          NGHPLS(2,IR,IP)=0
          NGHPLS(3,IR,IP)=0
          NGHPLS(4,IR,IP)=0
        ENDDO
      ENDDO

      DO IR=1,NR1STM
        DO K=1,NPPLG
          DO IP=NPOINT(1,K),NPOINT(2,K)
C  RADIAL NEIGHBORS
            NGHPOL(1,IR,IP)=IR-1
            NGHPOL(3,IR,IP)=MOD(IR+1,NR1ST)
            NGHPLS(1,IR,IP)=IR-1
            NGHPLS(3,IR,IP)=IR
C  POLOIDAL NEIGHBORS
            IF (IP.GT.NPOINT(1,K).AND.IP.LT.NPOINT(2,K)-1) THEN
C  INNER CELL
              NGHPOL(4,IR,IP)=IP-1
              NGHPOL(2,IR,IP)=IP+1
              NGHPLS(4,IR,IP)=IP-1
              NGHPLS(2,IR,IP)=IP
            ELSE IF (IP.EQ.NPOINT(1,K)) THEN
C  FIRST CELL OR SURFACE IN A PART OF A POLYGON
              NGHPOL(2,IR,IP)=IP+1
              NGHPLS(2,IR,IP)=IP
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
              DO IPART=1,NPPLG
                JP=NPOINT(2,IPART)
                IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
     .               (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS6 .AND.
     .              ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .               (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
                  NGHPOL(4,IR,IP)=JP-1
                  NGHPLS(4,IR,IP)=JP-1
                  GOTO 1
                ENDIF
              ENDDO
            ELSEIF (IP.EQ.NPOINT(2,K)-1) THEN
C  LAST CELL IN A PART OF A POLYGON
              NGHPOL(4,IR,IP)=IP-1
              NGHPLS(4,IR,IP)=IP-1
              NGHPLS(2,IR,IP)=IP
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
              DO IPART=1,NPPLG
                JP=NPOINT(1,IPART)
                IF (((XPOL(IR,IP+1)-XPOL(IR,JP))**2+
     .               (YPOL(IR,IP+1)-YPOL(IR,JP))**2).LT.EPS6 .AND.
     .              ((XPOL(IR+1,IP+1)-XPOL(IR+1,JP))**2+
     .               (YPOL(IR+1,IP+1)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
                  NGHPOL(2,IR,IP)=JP
                  GOTO 1
                ENDIF
              ENDDO
            ELSEIF (IP.EQ.NPOINT(2,K)) THEN
C  LAST SURFACE IN A PART OF A POLYGON
              NGHPLS(4,IR,IP)=IP-1
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
              DO IPART=1,NPPLG
                JP=NPOINT(1,IPART)
                IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
     .               (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS6 .AND.
     .              ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .               (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
                  NGHPLS(2,IR,IP)=JP
                  GOTO 1
                ENDIF
              ENDDO
            ENDIF
1         ENDDO
        ENDDO
      ENDDO

      baum => NewTree()
      NNODES=0
      DO IR=1,NR1STM
        DO K=1,NPPLG
          DO IP=NPOINT(1,K),NPOINT(2,K)-1
            IN=IR+(IP-1)*NR1ST
! IR, IP
            IC=NNODES+1
            inserted=.false.
            call insert (baum, xpol(ir,ip), ypol(ir,ip), 0._DP,
     .                   1._DP, ic, inserted)
            IF (INSERTED) THEN
              NNODES=NNODES+1
              XPOINT(NNODES)=XPOL(IR,IP)
              YPOINT(NNODES)=YPOL(IR,IP)
            END IF
            ALLOCATE (CUR)
            CUR%NOCELL = IN
            CUR%NEXT_CELL => COORCELL(IC)%PCELL
            COORCELL(IC)%PCELL => CUR
            INDPOINT(IR,IP) = IC
! IR+1, IP
            IC=NNODES+1
            inserted=.false.
            call insert (baum, xpol(ir+1,ip), ypol(ir+1,ip), 0._DP,
     .                   1._DP, ic, inserted)
            IF (INSERTED) THEN
              NNODES=NNODES+1
              XPOINT(NNODES)=XPOL(IR+1,IP)
              YPOINT(NNODES)=YPOL(IR+1,IP)
            END IF
            ALLOCATE (CUR)
            CUR%NOCELL = IN
            CUR%NEXT_CELL => COORCELL(IC)%PCELL
            COORCELL(IC)%PCELL => CUR
            INDPOINT(IR+1,IP) = IC
! IR+1, IP+1
            IC=NNODES+1
            inserted=.false.
            call insert (baum, xpol(ir+1,ip+1), ypol(ir+1,ip+1), 0._DP,
     .                   1._DP, ic, inserted)
            IF (INSERTED) THEN
              NNODES=NNODES+1
              XPOINT(NNODES)=XPOL(IR+1,IP+1)
              YPOINT(NNODES)=YPOL(IR+1,IP+1)
            END IF
            ALLOCATE (CUR)
            CUR%NOCELL = IN
            CUR%NEXT_CELL => COORCELL(IC)%PCELL
            COORCELL(IC)%PCELL => CUR
            INDPOINT(IR+1,IP+1) = IC
! IR, IP+1
            IC=NNODES+1
            inserted=.false.
            call insert (baum, xpol(ir,ip+1), ypol(ir,ip+1), 0._DP,
     .                   1._DP, ic, inserted)
            IF (INSERTED) THEN
              NNODES=NNODES+1
              XPOINT(NNODES)=XPOL(IR,IP+1)
              YPOINT(NNODES)=YPOL(IR,IP+1)
            END IF
            ALLOCATE (CUR)
            CUR%NOCELL = IN
            CUR%NEXT_CELL => COORCELL(IC)%PCELL
            COORCELL(IC)%PCELL => CUR
            INDPOINT(IR,IP+1) = IC
          ENDDO
        ENDDO
      ENDDO
      call DestroyTree(baum)

C     WRITE (iunout,*) '  IR    IP    S1    S2    S3    S4'
C     DO IR=1,NR1STM
C       DO IP=1,NP2NDM
C         WRITE (iunout,'(6I6)') IR,IP,(NGHPOL(K,IR,IP),K=1,4)
C         WRITE (iunout,'(6I6)') IR,IP,(NGHPLS(K,IR,IP),K=1,4)
C       ENDDO
C     ENDDO

      END




C ===== SOURCE: surtst.f
C
C
      SUBROUTINE SURTST(X,Y,Z,N,L)
C
C  THIS SUBROUTINE TESTS, WHETHER A POINT X,Y,Z FULLFILLS THE
C  THE BOUNDARY CONDITIONS FOR THE ADDITIONAL SURFACE NO. N OR NOT
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X, Y, Z
      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(OUT) :: L
      REAL(DP) :: XMS1, XMS2, XMS3, XLS1, XLS2, XLS3
      INTEGER :: I

      L=.TRUE.
      IF (RLB(N).EQ.1..OR.RLB(N).EQ.1.5) THEN
               L=L.AND.XLIMS1(1,N).LE.X.AND.X.LE.XLIMS2(1,N)
     .            .AND.YLIMS1(1,N).LE.Y.AND.Y.LE.YLIMS2(1,N)
     .            .AND.ZLIMS1(1,N).LE.Z.AND.Z.LE.ZLIMS2(1,N)
        IF (RLB(N).EQ.1.5) L=.NOT.L
      ELSEIF (RLB(N).LE.0.D0) THEN
        DO 1 I=1,ILIN(N)
           L=L.AND.
     .       ALIMS(I,N)+XLIMS(I,N)*X+YLIMS(I,N)*Y+ZLIMS(I,N)*Z.LE.0.
1       CONTINUE
        DO 2 I=1,ISCN(N)
           L=L.AND.
     .        ALIMS0(I,N)+XLIMS1(I,N)*X+YLIMS1(I,N)*Y+ZLIMS1(I,N)*Z
     .       +XLIMS2(I,N)*X*X+YLIMS2(I,N)*Y*Y+ZLIMS2(I,N)*Z*Z
     .       +XLIMS3(I,N)*X*Y+YLIMS3(I,N)*X*Z+ZLIMS3(I,N)*Y*Z.LE.0.
2       CONTINUE
      ELSEIF (RLB(N).GE.3.) THEN
        XMS1=X*PS13(1,N)+Y*PS13(2,N)+Z*PS13(3,N)+P1A(N)
        XLS1=X*PS23(1,N)+Y*PS23(2,N)+Z*PS23(3,N)+P2A(N)
        L=XMS1.GE.0..AND.XLS1.GE.0..AND.XMS1+XLS1.LE.1.
        IF (RLB(N).GE.4.AND..NOT.L) THEN
          XMS2=X*PS24(1,N)+Y*PS24(2,N)+Z*PS24(3,N)+P1B(N)
          XLS2=X*PS34(1,N)+Y*PS34(2,N)+Z*PS34(3,N)+P2B(N)
          L=XMS2.GE.0..AND.XLS2.GE.0..AND.XMS2+XLS2.LE.1.
          IF (RLB(N).GE.5.AND..NOT.L) THEN
            XMS3=X*PS35(1,N)+Y*PS35(2,N)+Z*PS35(3,N)+P1C(N)
            XLS3=X*PS45(1,N)+Y*PS45(2,N)+Z*PS45(3,N)+P2C(N)
            L=XMS3.GE.0..AND.XLS3.GE.0..AND.XMS3+XLS3.LE.1.
          ENDIF
        ENDIF
      ELSE
        WRITE (iunout,*) 'ERROR IN SUBROUTINE SURTST. RLB= ',RLB(N)
        CALL EXIT_OWN(1)
      ENDIF
      RETURN
      END
C ===== SOURCE: xshadd.f
C
C
      SUBROUTINE XSHADD (XSH,ILINI,ILEND)
C
C  SHIFT CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES IN X DIRECTION
C  OLD ORIGIN: X0=0.    (X,Y,Z)  SYSTEM
C  NEW ORIGIN: XN=-XSH   (X',Y,Z) SYSTEM
C
C     X'=X+XSH
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XSH
      INTEGER, INTENT(IN) :: ILINI, ILEND
      REAL(DP) :: XQ
      INTEGER :: I, J
C
      XQ=XSH*XSH
C
      DO 100 I=ILINI,ILEND
C
C
C   CHANGE COEFFICIENTS FOR ALGEBRAIC EQUATION
          A0LM(I)=A0LM(I)-A1LM(I)*XSH+A4LM(I)*XQ
          A1LM(I)=A1LM(I)-2.*A4LM(I)*XSH
          A2LM(I)=A2LM(I)-A7LM(I)*XSH
          A3LM(I)=A3LM(I)-A8LM(I)*XSH
          IF (RLB(I).LT.0.) THEN
C   CHANGE COEFFICIENTS IN LINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 10 J=1,ILIN(I)
              ALIMS(J,I)=ALIMS(J,I)-XLIMS(J,I)*XSH
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
              ALIMS0(J,I)=ALIMS0(J,I)-XLIMS1(J,I)*XSH+XLIMS2(J,I)*XQ
              XLIMS1(J,I)=XLIMS1(J,I)-2.*XLIMS2(J,I)*XSH
              YLIMS1(J,I)=YLIMS1(J,I)-XLIMS3(J,I)*XSH
              ZLIMS1(J,I)=ZLIMS1(J,I)-YLIMS3(J,I)*XSH
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER
            XLIMS1(1,I)=XLIMS1(1,I)+XSH
            XLIMS2(1,I)=XLIMS2(1,I)+XSH
C
C
          ELSE
C   POINT OPTIONS
            P1(1,I)=P1(1,I)+XSH
            P2(1,I)=P2(1,I)+XSH
            IF (RLB(I).GE.3.) P3(1,I)=P3(1,I)+XSH
            IF (RLB(I).GE.4.) P4(1,I)=P4(1,I)+XSH
            IF (RLB(I).GE.5.) P5(1,I)=P5(1,I)+XSH
            IF (RLB(I).GE.6.) P6(1,I)=P6(1,I)+XSH
          ENDIF
C
100   CONTINUE
C
      RETURN
      END
C ===== SOURCE: yshadd.f
C
C
      SUBROUTINE YSHADD (YSH,ILINI,ILEND)
C
C  SHIFT CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES IN Y DIRECTION
C  OLD ORIGIN: YO=O.    (X,Y,Z)  SYSTEM
C  NEW ORIGIN: YN=-YSH   (X,Y',Z) SYSTEM
C
C    Y'=Y+YSH
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: YSH
      INTEGER, INTENT(IN) :: ILINI, ILEND
      REAL(DP) :: YQ
      INTEGER :: I, J
      YQ=YSH*YSH
C
      DO 100 I=ILINI,ILEND
C
C
C   CHANGE COEFFICIENTS FOR ALGEBRAIC EQUATION
          A0LM(I)=A0LM(I)-A2LM(I)*YSH+A5LM(I)*YQ
          A2LM(I)=A2LM(I)-2.*A5LM(I)*YSH
          A1LM(I)=A1LM(I)-A7LM(I)*YSH
          A3LM(I)=A3LM(I)-A9LM(I)*YSH
          IF (RLB(I).LT.0.) THEN
C   CHANGE COEFFICIENTS IN LINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 10 J=1,ILIN(I)
              ALIMS(J,I)=ALIMS(J,I)-YLIMS(J,I)*YSH
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
              ALIMS0(J,I)=ALIMS0(J,I)-YLIMS1(J,I)*YSH+YLIMS2(J,I)*YQ
              YLIMS1(J,I)=YLIMS1(J,I)-2.*YLIMS2(J,I)*YSH
              XLIMS1(J,I)=XLIMS1(J,I)-XLIMS3(J,I)*YSH
              ZLIMS1(J,I)=ZLIMS1(J,I)-ZLIMS3(J,I)*YSH
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER
            YLIMS1(1,I)=YLIMS1(1,I)+YSH
            YLIMS2(1,I)=YLIMS2(1,I)+YSH
C
          ELSE
C
C   POINT OPTIONS
            P1(2,I)=P1(2,I)+YSH
            P2(2,I)=P2(2,I)+YSH
            IF (RLB(I).GE.3.) P3(2,I)=P3(2,I)+YSH
            IF (RLB(I).GE.4.) P4(2,I)=P4(2,I)+YSH
            IF (RLB(I).GE.5.) P5(2,I)=P5(2,I)+YSH
            IF (RLB(I).GE.6.) P6(2,I)=P6(2,I)+YSH
          ENDIF
C
100   CONTINUE
C
      RETURN
      END
C ===== SOURCE: zshadd.f
C
C
      SUBROUTINE ZSHADD (ZSH,ILINI,ILEND)
C
C  SHIFT CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES IN Z DIRECTION
C  OLD ORIGIN: ZO=O.    (X,Y,Z)  SYSTEM
C  NEW ORIGIN: ZN=-ZSH   (X,Y,Z') SYSTEM
C
C    Z'=Z+ZSH
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: ZSH
      INTEGER, INTENT(IN) :: ILINI, ILEND
      REAL(DP) :: ZQ
      INTEGER :: I, J

      ZQ=ZSH*ZSH
C
      DO 100 I=ILINI,ILEND
C
C
C   CHANGE COEFFICIENTS FOR ALGEBRAIC EQUATION
          A0LM(I)=A0LM(I)-A3LM(I)*ZSH+A6LM(I)*ZQ
          A3LM(I)=A3LM(I)-2.*A6LM(I)*ZSH
          A1LM(I)=A1LM(I)-A8LM(I)*ZSH
          A2LM(I)=A2LM(I)-A9LM(I)*ZSH
          IF (RLB(I).LT.0.) THEN
C   CHANGE COEFFICIENTS IN LINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 10 J=1,ILIN(I)
              ALIMS(J,I)=ALIMS(J,I)-ZLIMS(J,I)*ZSH
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
              ALIMS0(J,I)=ALIMS0(J,I)-ZLIMS1(J,I)*ZSH+ZLIMS2(J,I)*ZQ
              ZLIMS1(J,I)=ZLIMS1(J,I)-2.*ZLIMS2(J,I)*ZSH
              XLIMS1(J,I)=XLIMS1(J,I)-YLIMS3(J,I)*ZSH
              YLIMS1(J,I)=YLIMS1(J,I)-ZLIMS3(J,I)*ZSH
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER
            ZLIMS1(1,I)=ZLIMS1(1,I)+ZSH
            ZLIMS2(1,I)=ZLIMS2(1,I)+ZSH
C
          ELSE
C
C   POINT OPTIONS
            P1(3,I)=P1(3,I)+ZSH
            P2(3,I)=P2(3,I)+ZSH
            IF (RLB(I).GE.3.) P3(3,I)=P3(3,I)+ZSH
            IF (RLB(I).GE.4.) P4(3,I)=P4(3,I)+ZSH
            IF (RLB(I).GE.5.) P5(3,I)=P5(3,I)+ZSH
            IF (RLB(I).GE.6.) P6(3,I)=P6(3,I)+ZSH
          ENDIF
C
100   CONTINUE
C
      RETURN
      END
