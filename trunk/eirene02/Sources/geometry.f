c === ROUTINE: clltst
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
        WRITE (6,*) 'WRONG CELL-NUMBER IN RADIAL DIRECTION'
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
          WRITE (6,*) 'ERROR EXIT IN CLLTST, NLPOL ',NPCELL
          CALL EXIT
        ENDIF
        IF (NTEST1.NE.NPCELL.AND..NOT.NLSRFY) THEN
          WRITE (6,*) 'WRONG CELL-NUMBER IN POLOIDAL DIRECTION'
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
            WRITE (6,*) 'PHI,PHITEST  ',PHIT,PHITEST,PHIT-PHITEST,NTCELL
            CALL MASJ1 ('NPANU   ',NPANU)
          ENDIF
          NTEST2=LEARCA(PHIT,ZSURF,1,NTTRA,1,'CLLTST      ')
C
        ENDIF
        IF (NTEST2.NE.NTCELL.AND..NOT.NLSRFZ) THEN
          WRITE (6,*) 'WRONG CELL-NUMBER IN TOROIDAL DIRECTION'
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
          WRITE (6,*) 'PHI,PHITEST  ',PHIT,PHITEST,PHIT-PHITEST,NTEST2
          CALL MASJ1 ('NPANU   ',NPANU)
        ENDIF
C
        X0T=X0+RMTOR
        Z0T=Z0
        TTT=Z0T/(X0T*TANAL)
        IF (ABS(TTT).GT.1.+EPS10) THEN
          WRITE (6,*) 'WRONG CO-ORDINATES IN TOROIDAL DIRECTION'
          CALL MASJ1 ('NPANU   ',NPANU)
          CALL MASR3 ('X0,Z0,TTT                ',X0,Z0,TTT)
          RETURN 1
        ENDIF
C
        IF (NTEST2.NE.IPERID.AND..NOT.NLSRFZ) THEN
          WRITE (6,*) 'WRONG CELL-NUMBER IN TOROIDAL DIRECTION'
          CALL MASJ3 ('IPERID,NTEST2,NPANU     ',IPERID,NTEST2,NPANU)
        ENDIF
      ENDIF
C
      RETURN
      END
c === ROUTINE: fzrtor
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
            WRITE (6,*) 'ERROR IN FZRTOR: WRONG TOROIDAL CELL INDEX'
            WRITE (6,*) 'XOLD,ZOLD,XNEW,PH  ',XOLD,ZOLD,XNEW,PH
            WRITE (6,*) 'NNEW,NTST ',NNEW,NTEST
            WRITE (6,*) 'ABS(PH-ZSURF(NTST)),ABS(PH-ZSURF(NTST+1))'
            WRITE (6,*)  ABS(PH-ZSURF(NTEST)),ABS(PH-ZSURF(NTEST+1))
            WEIGHT=0.
            LGPART=.FALSE.
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
c === ROUTINE: fzrtra
C
C
      SUBROUTINE FZRTRA(X,Z,PH,NNEW)

      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGRID

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X, Z
      REAL(DP), INTENT(OUT) :: PH
      INTEGER, INTENT(OUT) :: NNEW
      REAL(DP) :: Z1, X01
      INTEGER :: LEARCA
C
C
C  FIND X,Z, NNEW,   FROM X,PH   (X=X??)
      Z1=ZSURF(1)
      PH=MOD(PH+PI2A-Z1,PI2A)+Z1
      NNEW=LEARCA(PH,ZSURF,1,NTTRA,1,'FZRTRA ')
      IF (NNEW.LE.0.OR.NNEW.GT.NTTRAM) THEN
        WRITE (6,*) 'NT OUT OF RANGE IN FZRTRA '
        WRITE (6,*) PH,ZHALF,NNEW
        CALL EXIT
      ENDIF
      X01=X+RMTOR
      CALL FZRTRI(X,Z,NNEW,X01,PH,NNEW)
      RETURN
      END
c === ROUTINE: frztri
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
c === ROUTINE: learc1
C
C
C      FUNCTION LEARC1 (X,Y,Z,IPO,IAN,IEN,LOGX,LOGY,NP,TEXT)
C      FUNCTION LEARC2 (X,Y,    IR,                ,NP,TEXT)
C      SUBROUTINE FZRTOR(X,Z,NOLD,XNEW,PH,NNEW,LTEST,NTEST)
C      SUBROUTINE FZRTRA(X,Z,PH,NNEW)
C      SUBROUTINE FZRTRI(X,Z,NNEW,XOLD,PH,NOLD)
C      FUNCTION AREAA (R,N,ARCA,YR,EP1R,ELLR)
C
C
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
C   IF .NOT.LOGX AND .NOT.LOGY
C     SEARCH IN RADIAL CELLS IAN AND IEN, I.E.
C     SEARCH BETWEEN (!!!) RADIAL SURFACES IAN AND IEN+1
C     THIS SEARCH COVERS THE WHOLE POLOIDAL RANGE
C   IF LOGX
C     SEARCH ON (!!!) RADIAL SURF. IAN FOR POLOIDAL MESH NUMBER IPO
C   IF LOGY
C     SEARCH ON (!!!) POLOIDAL SURF. IAN FOR RADIAL MESH NUMBER LEARC1
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
     .          D4, DELTAX, DELTAY
      REAL(DP), SAVE :: XMIN, YMIN, DISTX, DISTY, XMAX, YMAX, 
     .                  EPDY, EPDXDY, EPDX
      INTEGER, SAVE :: IFIRST
      INTEGER :: K, L, IM, LM, IEP, KH, LEARCT, LEAUSR, IMARK, LMARK,
     .           I, J, IE, LEARC1, IA,  
     .           IX, IY, INUM, IHEADX1, IHEADX2, IHEADY1, IHEADY2
c slmod begin
      INTEGER :: IR1, VP1
c slmod end
      REAL(DP), ALLOCATABLE, SAVE ::
     R D12(:,:), D12I(:,:), D14(:,:), D14I(:,:)
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
!pb      SAVE
      DATA IFIRST /0/
C
c slmod begin (debug)
c...  Does the broad search really not work?

      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(4X,A)')
     .    'LEARC1: (X,Y,IPO,IAN,IEN,LOGX,LOGY,NP,TEXT)'
        WRITE(6,'(10X,2F9.3,2X,I3,2X,2I3,2X,2L2,2X,I4,2X,A)')
     .    x,y,ipo,ian,ien,logx,logy,np,text

        IF (.NOT.LOGX.AND..NOT.LOGY) THEN
c          WRITE(6,'(4X,A)') 'LEARC1: *** ERROR! BROAD SEARCH DOES '//
c     .                      'NOT WORK ***'
c          WRITE(0,'(4X,A)') 'LEARC1: *** ERROR! BROAD SEARCH DOES '//
c     .                      'NOT WORK ***'
        ENDIF
      ENDIF

      IF (gridopt.EQ.1.AND.text.NE.'STDCOL'.AND.
     .                     text.NE.'ADDCOL  '.AND.
     .                     text.NE.'TORCOL 1') THEN
        WRITE(0,*) 'WARNING (Learc1): Function call of unknown ',
     .             'origin'
        WRITE(0,*) '     TEXT = ''',text,''''

        WRITE(6,*) 'WARNING (Learc1): Function call of unknown ',
     .             'origin'
        WRITE(6,*) '     TEXT = ''',text,''''
      ENDIF
c slmod end
      IA=IAN
      IE=IEN
C
      IF (LEVGEO.EQ.4) THEN
C
        IF (IFIRST.EQ.0) THEN
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
          XMIN = XMIN *(1.-EPS5)
          XMAX = XMAX *(1.+EPS5)
          YMIN = YMIN *(1.-EPS5)
          YMAX = YMAX *(1.+EPS5)
          DISTX=(XMAX-XMIN)/100.D0
          DISTY=(YMAX-YMIN)/100.D0
          EPDX=DISTX*EPS5
          EPDY=DISTY*EPS5
          EPDXDY=EPDX*EPDY
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
          DO 1 I=1,NR1ST
            DO 2 L=1,NP2NDM
c slmod begin
              IF (GRIDOPT.EQ.1) THEN
c...            The use of DIVSUR is a temporary fix.  There really needs to be some
c               parameter in the non-standard surface specifications that identifies
c               which radial surface of the ring the surface corresponds to:
                IF (I.GT.DIVSUR) THEN
                  IR1=I-1
                  VP1=2
                ELSE  
                  IR1=I
                  VP1=1
                ENDIF
                D12(I,L)=SQRT((XVERT(IR1,L,VP1)-XVERT(IR1,L+1,VP1))**2+
     .                        (YVERT(IR1,L,VP1)-YVERT(IR1,L+1,VP1))**2)
c...            Temp:
                DUM = SQRT((XPOL(I,L)-XPOL(I,L+1))**2+
     .                     (YPOL(I,L)-YPOL(I,L+1))**2)
                CALL CHECKNUM('zick1',i,l,d12(i,l),DUM)
              ELSE
                D12(I,L)=SQRT((XPOL(I,L)-XPOL(I,L+1))**2+
     .                        (YPOL(I,L)-YPOL(I,L+1))**2)
              ENDIF
c
c              D12(I,L)=SQRT((XPOL(I,L)-XPOL(I,L+1))**2+
c     .                      (YPOL(I,L)-YPOL(I,L+1))**2)
c slmod end
              D12I(I,L)=1./(ABS(D12(I,L))+EPS60)
2           CONTINUE
1         CONTINUE
          DO 3 I=1,NR1STM
            DO 4 L=1,NP2ND
c slmod begin
              IF (GRIDOPT.EQ.1) THEN
                D14(I,L)=SQRT((XVERT(I,L,1)-XVERT(I,L,2))**2+
     .                        (YVERT(I,L,1)-YVERT(I,L,2))**2)
c...            Temp:
                DUM = SQRT((XPOL(I,L)-XPOL(I+1,L))**2+
     .                     (YPOL(I,L)-YPOL(I+1,L))**2)
                CALL CHECKNUM('ZICK2',I,L,D14(I,L),DUM)
              ELSE
                D14(I,L)=SQRT((XPOL(I,L)-XPOL(I+1,L))**2+
     .                        (YPOL(I,L)-YPOL(I+1,L))**2)
              ENDIF
c
c              D14(I,L)=SQRT((XPOL(I,L)-XPOL(I+1,L))**2+
c     .                      (YPOL(I,L)-YPOL(I+1,L))**2)
c slmod end
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
c slmod begin
              IF (GRIDOPT.EQ.1) THEN
                IF (I.GT.DIVSUR) THEN
                  IR1=I-1
                  VP1=2
                ELSE  
                  IR1=I
                  VP1=1
                ENDIF
                XMIN = MIN(XMIN,XVERT(IR1,L,VP1))
                YMIN = MIN(YMIN,YVERT(IR1,L,VP1))
                XMAX = MAX(XMAX,XVERT(IR1,L,VP1))
                YMAX = MAX(YMAX,YVERT(IR1,L,VP1))
c...            Temp:
                DUM = MIN(XMIN,XPOL(I,L))
                CALL CHECKNUM('TICK1',I,L,XMIN,DUM)
                DUM = MIN(YMIN,YPOL(I,L))
                CALL CHECKNUM('TICK2',I,L,YMIN,DUM)
                DUM = MAX(XMAX,XPOL(I,L))
                CALL CHECKNUM('TICK3',I,L,XMAX,DUM)
                DUM = MAX(YMAX,YPOL(I,L))
                CALL CHECKNUM('TICK4',I,L,YMAX,DUM)
              ELSE
                XMIN = MIN(XMIN,XPOL(I,L))
                YMIN = MIN(YMIN,YPOL(I,L))
                XMAX = MAX(XMAX,XPOL(I,L))
                YMAX = MAX(YMAX,YPOL(I,L))
              ENDIF
c
c              XMIN = MIN(XMIN,XPOL(I,L))
c              YMIN = MIN(YMIN,YPOL(I,L))
c              XMAX = MAX(XMAX,XPOL(I,L))
c              YMAX = MAX(YMAX,YPOL(I,L))
c slmod end
            ENDDO
          ENDDO
          XMIN = XMIN *(1.-EPS5)
          XMAX = XMAX *(1.+EPS5)
          YMIN = YMIN *(1.-EPS5)
          YMAX = YMAX *(1.+EPS5)
          DISTX=(XMAX-XMIN)/100.
          DISTY=(YMAX-YMIN)/100.
!pb          EPDX=DISTX*EPS5
!pb          EPDY=DISTY*EPS5
!pb          EPDXDY=EPDX*EPDY
          EPDX=DISTX*EPS10
          EPDY=DISTY*EPS10
          EPDXDY=EPDX+EPDY
          write (6,*) ' epdx, epdy, epdxdy ',epdx, epdy, epdxdy
C  FOR EACH POLYGON CELL (I,L) FIND THE RANGE IHEADX1,....IHEADY2
C                              SUCH THAT THIS CELL (I,L) IS ENTIRELY
C                              IN THAT SECTION OF THE REGULAR IX,IY GRID
          DO I=1,NR1STM
            DO K=1,NPPLG
            DO L=NPOINT(1,K),NPOINT(2,K)-1
c slmod begin
              IF (GRIDOPT.EQ.1) THEN
                XCMIN=MIN(XVERT(I,L  ,1),XVERT(I,L  ,2),
     .                    XVERT(I,L+1,2),XVERT(I,L+1,1))
                XCMAX=MAX(XVERT(I,L  ,1),XVERT(I,L  ,2),
     .                    XVERT(I,L+1,2),XVERT(I,L+1,1))
                YCMIN=MIN(YVERT(I,L  ,1),YVERT(I,L  ,2),
     .                    YVERT(I,L+1,2),YVERT(I,L+1,1))
                YCMAX=MAX(YVERT(I,L  ,1),YVERT(I,L  ,2),
     .                    YVERT(I,L+1,2),YVERT(I,L+1,1))
c...            Temp:
                DUM = MIN(XPOL(I  ,L  ),XPOL(I+1,L  ),
     .                    XPOL(I+1,L+1),XPOL(I  ,L+1))
                CALL CHECKNUM('DICK1',I,L,XCMIN,DUM)
                DUM = MAX(XPOL(I  ,L  ),XPOL(I+1,L  ),
     .                    XPOL(I+1,L+1),XPOL(I  ,L+1))
                CALL CHECKNUM('DICK2',I,L,XCMAX,DUM)
                DUM = MIN(YPOL(I  ,L  ),YPOL(I+1,L  ),
     .                    YPOL(I+1,L+1),YPOL(I  ,L+1))
                CALL CHECKNUM('DICK3',I,L,YCMIN,DUM)
                DUM = MAX(YPOL(I  ,L  ),YPOL(I+1,L  ),
     .                    YPOL(I+1,L+1),YPOL(I  ,L+1))
                CALL CHECKNUM('DICK4',I,L,YCMAX,DUM)

              ELSE
                XCMIN=MIN(XPOL(I  ,L  ),XPOL(I+1,L  ),
     .                    XPOL(I+1,L+1),XPOL(I  ,L+1))
                XCMAX=MAX(XPOL(I  ,L  ),XPOL(I+1,L  ),
     .                    XPOL(I+1,L+1),XPOL(I  ,L+1))
                YCMIN=MIN(YPOL(I  ,L  ),YPOL(I+1,L  ),
     .                    YPOL(I+1,L+1),YPOL(I  ,L+1))
                YCMAX=MAX(YPOL(I  ,L  ),YPOL(I+1,L  ),
     .                    YPOL(I+1,L+1),YPOL(I  ,L+1))
              ENDIF
c
c              XCMIN=MIN(XPOL(I,L),XPOL(I+1,L),XPOL(I+1,L+1),XPOL(I,L+1))
c              XCMAX=MAX(XPOL(I,L),XPOL(I+1,L),XPOL(I+1,L+1),XPOL(I,L+1))
c              YCMIN=MIN(YPOL(I,L),YPOL(I+1,L),YPOL(I+1,L+1),YPOL(I,L+1))
c              YCMAX=MAX(YPOL(I,L),YPOL(I+1,L),YPOL(I+1,L+1),YPOL(I,L+1))
c slmod end
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
c slmod begin
            IF (GRIDOPT.EQ.1) THEN
              X1=XVERT(I,L  ,1)
              Y1=YVERT(I,L  ,1)
              X2=XVERT(I,L+1,1)
              Y2=YVERT(I,L+1,1)
              X3=XVERT(I,L+1,2)
              Y3=YVERT(I,L+1,2)
c...          Temp:
              DUM = XPOL(I,L)
              CALL CHECKNUM('RICK1',I,L,X1,DUM)
              DUM = YPOL(I,L)
              CALL CHECKNUM('RICK2',I,L,Y1,DUM)
              DUM = XPOL(I,L+1)
              CALL CHECKNUM('RICK3',I,L,X2,DUM)
              DUM = YPOL(I,L+1)
              CALL CHECKNUM('RICK4',I,L,Y2,DUM)
              DUM = XPOL(I+1,L+1)
              CALL CHECKNUM('RICK5',I,L,X3,DUM)
              DUM = YPOL(I+1,L+1)
              CALL CHECKNUM('RICK6',I,L,Y3,DUM)
            ELSE
              X1=XPOL(I,L)
              Y1=YPOL(I,L)
              X2=XPOL(I,L+1)
              Y2=YPOL(I,L+1)
              X3=XPOL(I+1,L+1)
              Y3=YPOL(I+1,L+1)
            ENDIF
c
c            X1=XPOL(I,L)
c            Y1=YPOL(I,L)
c            X2=XPOL(I,L+1)
c            Y2=YPOL(I,L+1)
c            X3=XPOL(I+1,L+1)
c            Y3=YPOL(I+1,L+1)
c slmod end
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
c slmod begin
            IF (GRIDOPT.EQ.1) THEN
              X4=XVERT(I,L,2)
              Y4=YVERT(I,L,2)        
c...          Temp:
              DUM=XPOL(I+1,L)
              CALL CHECKNUM('RICK1',I,L,X4,DUM)
              DUM=YPOL(I+1,L)
              CALL CHECKNUM('RICK2',I,L,Y4,DUM)
            ELSE
              X4=XPOL(I+1,L)
              Y4=YPOL(I+1,L)
            ENDIF
c
c            X4=XPOL(I+1,L)
c            Y4=YPOL(I+1,L)
c slmod end
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
c slmod begin (debug)
            IF (printopt.EQ.2) 
     .        WRITE(6,'(4X,A,5I4)')
     .          'LEARC1: (I,IA,IEP,L,NP2NDM) ',i,ia,iep,l,np2ndm
c slmod end
121         IF (L.GT.NP2NDM) GOTO 21
c slmod begin
            IF (GRIDOPT.EQ.1) THEN
              IF (I.GT.1) THEN
                IR1=I-1
                VP1=2
              ELSE  
                IR1=I
                VP1=1
              ENDIF
c...          Do not process virtual cells:
              IF (RVRTAG(IR1,L).EQ.1) THEN
                ERR1=1.0D+30
                GOTO 21
              ENDIF
              XMX1=X-XVERT(IR1,L  ,VP1)
              YMY1=Y-YVERT(IR1,L  ,VP1)
              XMX2=X-XVERT(IR1,L+1,VP1)
              YMY2=Y-YVERT(IR1,L+1,VP1)
c...          Temp:
              CALL CHECKNUM('zing1',i,l  ,xmx1,X-XPOL(I,L  ))
              CALL CHECKNUM('zing2',i,l+1,xmx2,X-XPOL(I,L+1))
              CALL CHECKNUM('zing3',i,l  ,ymy1,Y-YPOL(I,L  ))
              CALL CHECKNUM('zing4',i,l  ,ymy2,Y-YPOL(I,L+1))
            ELSE
              XMX1=X-XPOL(I,L)
              YMY1=Y-YPOL(I,L)
              XMX2=X-XPOL(I,L+1)
              YMY2=Y-YPOL(I,L+1)
            ENDIF
c
c            XMX1=X-XPOL(I,L)
c            YMY1=Y-YPOL(I,L)
c            XMX2=X-XPOL(I,L+1)
c            YMY2=Y-YPOL(I,L+1)
c slmod end
            D1=SQRT(XMX1*XMX1+YMY1*YMY1)
            D2=SQRT(XMX2*XMX2+YMY2*YMY2)
            ERR1=ABS(D1+D2-D12(I,L))*D12I(I,L)
c slmod begin (debug)
            IF (printopt.EQ.2)
     .        WRITE(6,'(4X,A,4I5,1P,E15.8)')
     .          'LEARC1: Searching (I,L,IR1,VP1,ERR1) ',
     .          i,l,ir1,vp1,err1
c slmod end
            IF (ERR1.LT.ERRMIN) THEN
              IMARK=I
              LMARK=L
              ERRMIN=ERR1
            ENDIF
21          CONTINUE
c slmod begin
c...        Search both sides of each ring, after ring DIVSUR.  Again,
c           this is a temporary measure (see above):
            IF (I.GT.DIVSUR.AND.I.LE.NR1STM) THEN
c
c            IF (I .EQ. NR1STM) THEN
c slmod end
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
c slmod begin (debug)
        IF (printopt.EQ.2) WRITE(6,'(4X,A)') 'LEARC1: Checkpoint C'
c slmod end
        DO J=1,4
          HELPP => HELPCUR4(J)%P
          DO WHILE (ASSOCIATED(HELPP))
            I = HELPP%IX
            L = HELPP%IY
221         IF ((I .LT. IAN) .OR. (I .GT. IEN)) GOTO 22
c slmod begin
            IF (GRIDOPT.EQ.1) THEN
              XMX1=X-XVERT(I,L,1)
              YMY1=Y-YVERT(I,L,1)
              XMX4=X-XVERT(I,L,2)
              YMY4=Y-YVERT(I,L,2)
c...          Temp:
              CALL CHECKNUM('zing5',i,l,xmx1,X-XPOL(I  ,L))
              CALL CHECKNUM('zing6',i,l,xmx4,X-XPOL(I+1,L))
              CALL CHECKNUM('zing7',i,l,ymy1,Y-YPOL(I  ,L))
              CALL CHECKNUM('zing8',i,l,ymy4,Y-YPOL(I+1,L))
            ELSE
              XMX1=X-XPOL(I,L)
              YMY1=Y-YPOL(I,L)
              XMX4=X-XPOL(I+1,L)
              YMY4=Y-YPOL(I+1,L)
            ENDIF
c
c            XMX1=X-XPOL(I,L)
c            YMY1=Y-YPOL(I,L)
c            XMX4=X-XPOL(I+1,L)
c            YMY4=Y-YPOL(I+1,L)
c slmod end
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
c slmod begin
            IF (GRIDOPT.EQ.1) THEN
              XMX1=X-XVERT(I,IA,1)
              YMY1=Y-YVERT(I,IA,1)
              XMX4=X-XVERT(I,IA,2)
              YMY4=Y-YVERT(I,IA,2)
c...          Temp:
              CALL CHECKNUM('zing9 ',i,ia,xmx1,X-XPOL(I  ,IA))
              CALL CHECKNUM('zing10',i,ia,xmx4,X-XPOL(I+1,IA))
              CALL CHECKNUM('zing11',i,ia,ymy1,Y-YPOL(I  ,IA))
              CALL CHECKNUM('zing12',i,ia,ymy4,Y-YPOL(I+1,IA))
            ELSE
              XMX1=X-XPOL(I,IA)
              YMY1=Y-YPOL(I,IA)
              XMX4=X-XPOL(I+1,IA)
              YMY4=Y-YPOL(I+1,IA)
            ENDIF
c
c            XMX1=X-XPOL(I,IA)
c            YMY1=Y-YPOL(I,IA)
c            XMX4=X-XPOL(I+1,IA)
c            YMY4=Y-YPOL(I+1,IA)
c slmod end
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
        ELSEIF (INUM.GT.1.AND.ERRMIN.GT.EPS10) THEN
          CALL MASAGE ('WARNING FROM LEARC1, INUM.GT.1               ')
          CALL MASR2('X,Y             ',X,Y)
          WRITE (6,*) 'LEARC1 CALLED FROM SUBR. ',TEXT
          WRITE (6,*) 'ERRMIN= ',ERRMIN
          WRITE (6,*) 'NPANU,INUM,IM,LM= ',NP,INUM,IM,LM
          WRITE (6,*) 'IAN,IEN,LOGX,LOGY ',IAN,IEN,LOGX,LOGY
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
        CALL EXIT
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
        CALL EXIT
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
c === ROUTINE: learc2
C
C
      FUNCTION LEARC2(X,Y,NR,NP,TEXT)
C
C  THIS SUBROUTINE FINDS THE POLYGON INDEX "IPOLG"
C  ASSUMING THAT X,Y IS IN THE RADIAL ZONE NR
C
      USE PRECISION
      USE PARMMOD
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
c slmod begin (debug)
      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(4X,A                 )') 'LEARC2: (X,Y,NR,NP,TEXT)'
        WRITE(6,'(10X,2F9.3,2X,2I5,2X,A)')           x,y,nr,np,text
      ENDIF
c slmod end
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
c slmod begin
              IF (GRIDOPT.EQ.1) THEN
c...            Do not process virtual (zero volume) cells:
                IF (RVRTAG(NR,L).EQ.1) CYCLE
                X1(N,L)=XVERT(N,L  ,1)
                Y1(N,L)=YVERT(N,L  ,1)
                X2(N,L)=XVERT(N,L+1,1)
                Y2(N,L)=YVERT(N,L+1,1)
                X3(N,L)=XVERT(N,L+1,2)
                Y3(N,L)=YVERT(N,L+1,2)
                X4(N,L)=XVERT(N,L  ,2)
                Y4(N,L)=YVERT(N,L  ,2)
c...            Temp:
                CALL CHECKNUM('ting1',n,l,X1(N,L),XPOL(N,L)  )
                CALL CHECKNUM('ting2',n,l,Y1(N,L),YPOL(N,L)  )
                CALL CHECKNUM('ting3',n,l,X2(N,L),XPOL(N,L+1))
                CALL CHECKNUM('ting4',n,l,Y2(N,L),YPOL(N,L+1))
                CALL CHECKNUM('ting5',n,l,X3(N,L),XPOL(N+1,L+1))
                CALL CHECKNUM('ting6',n,l,Y3(N,L),YPOL(N+1,L+1))
                CALL CHECKNUM('ting7',n,l,X4(N,L),XPOL(N+1,L))
                CALL CHECKNUM('ting8',n,l,Y4(N,L),YPOL(N+1,L))
              ELSE
                X1(N,L)=XPOL(N,L)
                Y1(N,L)=YPOL(N,L)
                X2(N,L)=XPOL(N,L+1)
                Y2(N,L)=YPOL(N,L+1)
                X3(N,L)=XPOL(N+1,L+1)
                Y3(N,L)=YPOL(N+1,L+1)
                X4(N,L)=XPOL(N+1,L)
                Y4(N,L)=YPOL(N+1,L)
              ENDIF
c
c              X1(N,L)=XPOL(N,L)
c              Y1(N,L)=YPOL(N,L)
c              X2(N,L)=XPOL(N,L+1)
c              Y2(N,L)=YPOL(N,L+1)
c              X3(N,L)=XPOL(N+1,L+1)
c              Y3(N,L)=YPOL(N+1,L+1)
c              X4(N,L)=XPOL(N+1,L)
c              Y4(N,L)=YPOL(N+1,L)
c slmod end
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
c slmod begin
c...    Do not process virtual (zero volume) cells:
        IF (GRIDOPT.EQ.1.AND.RVRTAG(NR,L).EQ.1) CYCLE
c slmod end
        IF (ERR4(L).GT.1.D30) THEN
c slmod begin
          IF (GRIDOPT.EQ.1) THEN
            X2N=XVERT(NR,L  ,1)
            Y2N=YVERT(NR,L  ,1)
            X3N=XVERT(NR,L+1,1)
            Y3N=YVERT(NR,L+1,1)
            X4N=XVERT(NR,L+1,2)
            Y4N=YVERT(NR,L+1,2)
            X1N=XVERT(NR,L  ,2)
            Y1N=YVERT(NR,L  ,2)
c...        Temp:
            CALL CHECKNUM('jing1',nr,l,X2N,XPOL(NR,L))  
            CALL CHECKNUM('jing2',nr,l,Y2N,YPOL(NR,L))  
            CALL CHECKNUM('jing3',nr,l,X3N,XPOL(NR,L+1))
            CALL CHECKNUM('jing4',nr,l,Y3N,YPOL(NR,L+1))
            CALL CHECKNUM('jing5',nr,l,X4N,XPOL(NR+1,L+1))
            CALL CHECKNUM('jing6',nr,l,Y4N,YPOL(NR+1,L+1))
            CALL CHECKNUM('jing7',nr,l,X1N,XPOL(NR+1,L))
            CALL CHECKNUM('jing8',nr,l,Y1N,YPOL(NR+1,L))
          ELSE
            X2N=XPOL(NR,L)
            Y2N=YPOL(NR,L)
            X3N=XPOL(NR,L+1)
            Y3N=YPOL(NR,L+1)
            X4N=XPOL(NR+1,L+1)
            Y4N=YPOL(NR+1,L+1)
            X1N=XPOL(NR+1,L)
            Y1N=YPOL(NR+1,L)
          ENDIF
c
c          X2N=XPOL(NR,L)
c          Y2N=YPOL(NR,L)
c          X3N=XPOL(NR,L+1)
c          Y3N=YPOL(NR,L+1)
c          X4N=XPOL(NR+1,L+1)
c          Y4N=YPOL(NR+1,L+1)
c          X1N=XPOL(NR+1,L)
c          Y1N=YPOL(NR+1,L)
c slmod end
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
          ERR1(L)=ERRMIN
          ERR3(L)=ERRMIN
        ENDIF
        DO 110 L=NPOINT(1,K),NPOINT(2,K)-1
c slmod begin
c...      Do not process virtual (zero volume) cells:
          IF (GRIDOPT.EQ.1.AND.RVRTAG(NR,L).EQ.1) CYCLE
c slmod end
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
        WRITE (6,*) 'LEARC2 CALLED FROM SUBR. ',TEXT
        WRITE (6,*) 'ERRMIN= ',ERRMIN
        WRITE (6,*) 'NPANU,IM,LM= ',NP,IM,LM
      ELSEIF (INUM.GT.1.AND.ERRMIN.GT.EPS5) THEN
        CALL MASAGE ('WARNING FROM LEARC2, INUM.GT.1               ')
        CALL MASR2('X,Y             ',X,Y)
        WRITE (6,*) 'LEARC2 CALLED FROM SUBR. ',TEXT
        WRITE (6,*) 'ERRMIN= ',ERRMIN
        WRITE (6,*) 'NPANU,INUM,IM,LM= ',NP,INUM,IM,LM
      ENDIF
C
      LEARC2=LM
C
      RETURN
      END
c === ROUTINE: rotadd
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
C           WRITE (6,*) 'INVARIANTS FROM ROTADD: I= ',I
C           WRITE (6,*) 'BEFORE: S,DEL,T= ',S,DEL,T
C           WRITE (6,*) ' A0...A9 ',A0LM(I),A1LM(I),A2LM(I),A3LM(I),
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
C         WRITE (6,*) 'AFTER:  S,DEL,T= ',S,DEL,T
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
              WRITE (6,*) 'ROTADD: TO BE WRITTEN, RLB=1.5'
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
c === ROUTINE: sneigh
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

      IMPLICIT NONE

      INTEGER :: IR, IP, IPART, JP, K
c slmod begin
      REAL(DP) :: DIST1, DIST2
c slmod end

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
c slmod begin
                IF (GRIDOPT.EQ.1) THEN
                  DIST1=(XVERT(IR,IP,1)-XVERT(IR,JP,1))**2+
     .                  (YVERT(IR,IP,1)-YVERT(IR,JP,1))**2
                  DIST2=(XVERT(IR,IP,2)-XVERT(IR,JP,2))**2+
     .                  (YVERT(IR,IP,2)-YVERT(IR,JP,2))**2
c...              Temp:
                  DUM=(XPOL(IR,IP)-XPOL(IR,JP))**2+
     .                (YPOL(IR,IP)-YPOL(IR,JP))**2
                  CALL CHECKNUM('DAKA1',IR,IP,DIST1,DUM)
                  DUM=(XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .                (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2
                  CALL CHECKNUM('DAKA1',IR,IP,DIST2,DUM)

                  IF (DIST1.LT.EPS6.AND.DIST2.LT.EPS6) THEN
                    NGHPOL(4,IR,IP)=JP-1
                    NGHPLS(4,IR,IP)=JP-1
                    GOTO 1
                  ENDIF
                ELSE
                  IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
     .                 (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS6 .AND.
     .                ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .                 (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
                    NGHPOL(4,IR,IP)=JP-1
                    NGHPLS(4,IR,IP)=JP-1
                    GOTO 1
                  ENDIF
                ENDIF
c
c                IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
c     .               (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS6 .AND.
c     .              ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
c     .               (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
c                  NGHPOL(4,IR,IP)=JP-1
c                  NGHPLS(4,IR,IP)=JP-1
c                  GOTO 1
c                ENDIF
c slmod end
              ENDDO
            ELSEIF (IP.EQ.NPOINT(2,K)-1) THEN
C  LAST CELL IN A PART OF A POLYGON
              NGHPOL(4,IR,IP)=IP-1
              NGHPLS(4,IR,IP)=IP-1
              NGHPLS(2,IR,IP)=IP
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
              DO IPART=1,NPPLG
                JP=NPOINT(1,IPART)
c slmod begin
                IF (GRIDOPT.EQ.1) THEN
                  DIST1=(XVERT(IR,IP+1,1)-XVERT(IR,JP,1))**2+
     .                  (YVERT(IR,IP+1,1)-YVERT(IR,JP,1))**2
                  DIST2=(XVERT(IR,IP+1,2)-XVERT(IR,JP,2))**2+
     .                  (YVERT(IR,IP+1,2)-YVERT(IR,JP,2))**2
c...              Temp:
                  DUM=(XPOL(IR,IP+1)-XPOL(IR,JP))**2+
     .                (YPOL(IR,IP+1)-YPOL(IR,JP))**2
                  CALL CHECKNUM('JAKA1',IR,IP,DIST1,DUM)
                  DUM=(XPOL(IR+1,IP+1)-XPOL(IR+1,JP))**2+
     .                (YPOL(IR+1,IP+1)-YPOL(IR+1,JP))**2
                  CALL CHECKNUM('JAKA1',IR,IP,DIST2,DUM)

                  IF (DIST1.LT.EPS6.AND.DIST2.LT.EPS6) THEN
                    NGHPOL(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ELSE
                  IF   (((XPOL(IR,IP+1)-XPOL(IR,JP))**2+
     .                 (YPOL(IR,IP+1)-YPOL(IR,JP))**2).LT.EPS6 .AND.
     .                ((XPOL(IR+1,IP+1)-XPOL(IR+1,JP))**2+
     .                 (YPOL(IR+1,IP+1)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
                    NGHPOL(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ENDIF
c
c                  IF   (((XPOL(IR,IP+1)-XPOL(IR,JP))**2+
c     .                 (YPOL(IR,IP+1)-YPOL(IR,JP))**2).LT.EPS6 .AND.
c     .                ((XPOL(IR+1,IP+1)-XPOL(IR+1,JP))**2+
c     .                 (YPOL(IR+1,IP+1)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
c                    NGHPOL(2,IR,IP)=JP
c                    GOTO 1
c                  ENDIF
c slmod end
              ENDDO
            ELSEIF (IP.EQ.NPOINT(2,K)) THEN
C  LAST SURFACE IN A PART OF A POLYGON
              NGHPLS(4,IR,IP)=IP-1
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
              DO IPART=1,NPPLG
                JP=NPOINT(1,IPART)
c slmod begin
                IF (GRIDOPT.EQ.1) THEN
                  DIST1=(XVERT(IR,IP,1)-XVERT(IR,JP,1))**2+
     .                  (YVERT(IR,IP,1)-YVERT(IR,JP,1))**2
                  DIST2=(XVERT(IR,IP,2)-XVERT(IR,JP,2))**2+
     .                  (YVERT(IR,IP,2)-YVERT(IR,JP,2))**2
c...              Temp:
                  DUM=(XPOL(IR,IP)-XPOL(IR,JP))**2+
     .                (YPOL(IR,IP)-YPOL(IR,JP))**2
                  CALL CHECKNUM('MAKA1',IR,IP,DIST1,DUM)
                  DUM=(XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .                (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2
                  CALL CHECKNUM('MAKA1',IR,IP,DIST2,DUM)

                  IF (DIST1.LT.EPS6.AND.DIST2.LT.EPS6) THEN
                    NGHPLS(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ELSE
                  IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
     .                 (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS6 .AND.
     .                ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .                 (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
                    NGHPLS(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ENDIF
c
c                IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
c     .               (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS6 .AND.
c     .              ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
c     .               (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS6) THEN
c                  NGHPLS(2,IR,IP)=JP
c                  GOTO 1
c                ENDIF
c slmod end
              ENDDO
            ENDIF
1         ENDDO
        ENDDO
      ENDDO

C     WRITE (6,*) '  IR    IP    S1    S2    S3    S4'
C     DO IR=1,NR1STM
C       DO IP=1,NP2NDM
C         WRITE (6,'(6I6)') IR,IP,(NGHPOL(K,IR,IP),K=1,4)
C         WRITE (6,'(6I6)') IR,IP,(NGHPLS(K,IR,IP),K=1,4)
C       ENDDO
C     ENDDO
c slmod begin (sl)
      WRITE(0,*) 'CHECK CONNECTION MAP AND BE SURE THAT RADIAL '//
     .           'TRANSPORT IS OKAY'
c slmod end
      END
c === ROUTINE: surtst
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
        WRITE (6,*) 'ERROR IN SUBROUTINE SURTST. RLB= ',RLB(N)
        CALL EXIT
      ENDIF
      RETURN
      END
c === ROUTINE: timea
C
C
      SUBROUTINE TIMEA
C
C   1 ST INTERSECTION OF THE RAY X+T*VX,Y+T*VY,Z+T*VZ WITH ONE OF THE
C   ADDITIONAL SURFACES, DEFINED BY 2.ND ORDER EQUATIONS
C   IT IS ALSO CHECKED, WETHER THIS INTERSECTION TAKES PLACE INSIDE THE
C   SPECIFIED BOUNDARIES OF THOSE SURFACES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CLGIN
      USE CTRIG

      IMPLICIT NONE
C
      REAL(DP) :: XB(3), XC(3), XD(3), XE(3), XF(3), XG(3)
      REAL(DP) :: TADD, TMIN, XXR, V, VSAVE, A1, A2, PPP, ROT, VX, VY,
     .          VZ, VXS, VYS, VZS,VVS, TS, T, X, Y, Z, A3, VXJ, VYJ,
     .          VZJ, XLS2, XLS3, XMS3, XMS2, XX, YY, ZZ, XR, YR, ZR,
     .          XN, YN, ZN, TMA, TMI, WR, TST, XMS1, XLS1, TMX, TUP,
     .          DZ3, TL, XS, DZ, TEST4, TEST5, DY, DL, DX1, DX, SG,
     .          VV, VXX, VYY, VZZ, TMT, B, XK5, V445, XN5, V335, XK4,
     .          V334, V3V45, XN4, YS, ZS, AT, XNORM, V224, CG1, CG2,
     .          CG3, AR1, C, D, E, F, G, XK3, HELP, V2V34, XN3, V223
     .          AR3, AR2, V113, V1V23, V223, AR3, DY2
      INTEGER :: NN, NLLLI, NTNEW, NNTCLS, NNR, ICOUNT, IPLGN, ITRII,
     .           ISTS, I1000, J, IDEZ, I, IPERID, NTCELL, NCELL, MSURF,
     .           NLI, NLE, MASURF, ISPZ, NNTCL, LM2, LM1, IAB, JUM
c slmod begin
      INTEGER :: IREG
c slmod end
      LOGICAL :: LGJ, NLTRC, BITGET
      SAVE
C
      ENTRY TIMEA0
C
      IF (NLIMI.LT.1) RETURN
C
C
      DO 2 J=1,NLIMI
        IF (IGJUM0(J).NE.0) THEN
          IF (NLIMPB >= NLIMPS) THEN
            DO 1 I=0,NLIMI
              IGJUM1(I,J)=1
1           CONTINUE
          ELSE
            DO I=0,NLIMI
              CALL BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
            END DO
          END IF
        ENDIF
C
        ISWICH(1,J)=IDEZ(ILSWCH(J),1,6)
        IF (ISWICH(1,J).EQ.1) ISWICH(1,J)=-1
        IF (ISWICH(1,J).EQ.2) ISWICH(1,J)=1
        ISWICH(2,J)=IDEZ(ILSWCH(J),2,6)
        IF (ISWICH(2,J).EQ.1) ISWICH(2,J)=-1
        IF (ISWICH(2,J).EQ.2) ISWICH(2,J)=1
        ISWICH(3,J)=IDEZ(ILSWCH(J),3,6)
        IF (ISWICH(3,J).EQ.1) ISWICH(3,J)=-1
        IF (ISWICH(3,J).EQ.2) ISWICH(3,J)=1
        ISWICH(4,J)=IDEZ(ILSWCH(J),4,6)
        IF (ISWICH(4,J).EQ.1) ISWICH(4,J)=-1
        IF (ISWICH(4,J).EQ.2) ISWICH(4,J)=1
        ISWICH(5,J)=IDEZ(ILSWCH(J),5,6)
        IF (ISWICH(5,J).EQ.1) ISWICH(5,J)=-1
        IF (ISWICH(5,J).EQ.2) ISWICH(5,J)=1
        ISWICH(6,J)=IDEZ(ILSWCH(J),6,6)
        IF (ISWICH(6,J).EQ.1) ISWICH(6,J)=-1
        IF (ISWICH(6,J).EQ.2) ISWICH(6,J)=1
C
        IF (ISWICH(4,J).NE.0.OR.ISWICH(5,J).NE.0.OR.ISWICH(6,J).NE.0)
     .  THEN
          ILBLCK(J)=IDEZ(ILCELL(J),4,4)
          I1000=1000*ILBLCK(J)
          ILACLL(J)=ILCELL(J)-I1000
        ENDIF
2     CONTINUE
C
      DO 97 J=1,NLIMI
        IF (IGJUM0(J).NE.0) THEN
          IF (LEVGEO.EQ.4) THEN
            DO ITRII=1,NTRII
            DO IPLGN=1,3
              ISTS=ABS(INMTI(IPLGN,ITRII))
              IF (J.EQ.ISTS) GOTO 85
            ENDDO
            ENDDO
          ENDIF
          GOTO 97
        ENDIF
85      IF (RLB(J).LT.2.0) THEN
C
C   SURFACE COEFFICIENTS ARE INPUT
C
C  INDICATE INDEPENDENCE OF SURFACE EQUATION FROM X, Y, Z, RESP.
          IF (A1LM(J).EQ.0..AND.A4LM(J).EQ.0..AND.
     .        A7LM(J).EQ.0..AND.A8LM(J).EQ.0.)
     .    P3(1,J)=1.D55
          IF (A2LM(J).EQ.0..AND.A5LM(J).EQ.0..AND.
     .        A7LM(J).EQ.0..AND.A9LM(J).EQ.0.)
     .    P3(2,J)=1.D55
          IF (A3LM(J).EQ.0..AND.A6LM(J).EQ.0..AND.
     .        A8LM(J).EQ.0..AND.A9LM(J).EQ.0.)
     .    P3(3,J)=1.D55
          GOTO 90
        ELSEIF (RLB(J).GE.2) THEN
C
C  PLANE SURFACE, VERTICES ARE INPUT, SET SURFACE COEFFICIENTS:
C
          A4LM(J)=0.
          A5LM(J)=0.
          A6LM(J)=0.
          A7LM(J)=0.
          A8LM(J)=0.
          A9LM(J)=0.
        ENDIF
C
        IF (RLB(J).GE.2.6.AND.RLB(J).LT.2.9) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE X CO-ORDINATE
          IF (RLB(J).LT.2.75) RLB(J)=1.
          IF (RLB(J).GE.2.75) RLB(J)=1.5
          XLIMS1(1,J)=P1(1,J)
          XLIMS2(1,J)=P2(1,J)
          DZ3=(P1(3,J)-P2(3,J))
          DY2=(P1(2,J)-P2(2,J))
          A0LM(J)=DZ3*P1(2,J)-DY2*P1(3,J)
          A1LM(J)=0.
          A2LM(J)=-DZ3
          A3LM(J)=DY2
          YLIMS1(1,J)=MIN(P1(2,J),P2(2,J))
          YLIMS2(1,J)=MAX(P1(2,J),P2(2,J))
          ZLIMS1(1,J)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(1,J)=MAX(P1(3,J),P2(3,J))
          IF (P1(2,J).EQ.P2(2,J)) THEN
            YLIMS1(1,J)=YLIMS1(1,J)-0.1
            YLIMS2(1,J)=YLIMS2(1,J)+0.1
          ENDIF
          IF (P1(3,J).EQ.P2(3,J)) THEN
            ZLIMS1(1,J)=ZLIMS1(1,J)-0.1
            ZLIMS2(1,J)=ZLIMS2(1,J)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND X-INDEPENDENCE OF SURFACE EQUATION
          P3(1,J)=1.D55
          DL=SQRT(DY2**2+DZ3**2)
          DX=XLIMS2(1,J)-XLIMS1(1,J)
          IF (DX.GT.1.D20) DX=XDF
          SAREA(J)=DL*DX
          GOTO 90
        ELSEIF (RLB(J).GE.2.3.AND.RLB(J).LT.2.6) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE Y CO-ORDINATE
          IF (RLB(J).LT.2.45) RLB(J)=1.
          IF (RLB(J).GE.2.45) RLB(J)=1.5
          YLIMS1(1,J)=P1(2,J)
          YLIMS2(1,J)=P2(2,J)
          DZ3=(P1(3,J)-P2(3,J))
          DX1=(P1(1,J)-P2(1,J))
          A0LM(J)=DZ3*P1(1,J)-DX1*P1(3,J)
          A1LM(J)=-DZ3
          A2LM(J)=0.
          A3LM(J)=DX1
          XLIMS1(1,J)=MIN(P1(1,J),P2(1,J))
          XLIMS2(1,J)=MAX(P1(1,J),P2(1,J))
          ZLIMS1(1,J)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(1,J)=MAX(P1(3,J),P2(3,J))
          IF (P1(1,J).EQ.P2(1,J)) THEN
            XLIMS1(1,J)=XLIMS1(1,J)-0.1
            XLIMS2(1,J)=XLIMS2(1,J)+0.1
          ENDIF
          IF (P1(3,J).EQ.P2(3,J)) THEN
            ZLIMS1(1,J)=ZLIMS1(1,J)-0.1
            ZLIMS2(1,J)=ZLIMS2(1,J)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND Y-INDEPENDENCE OF SURFACE EQUATION
          P3(2,J)=1.D55
          DL=SQRT(DX1**2+DZ3**2)
          DY=YLIMS2(1,J)-YLIMS1(1,J)
          IF (DY.GT.1.E20) DY=YDF
          SAREA(J)=DL*DY
          GOTO 90
        ELSEIF (RLB(J).GE.2.0.AND.RLB(J).LT.2.3) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE Z CO-ORDINATE
          IF (RLB(J).LT.2.15) RLB(J)=1.
          IF (RLB(J).GE.2.15) RLB(J)=1.5
          ZLIMS1(1,J)=P1(3,J)
          ZLIMS2(1,J)=P2(3,J)
          DY2=(P1(2,J)-P2(2,J))
          DX1=(P1(1,J)-P2(1,J))
          A0LM(J)=DY2*P1(1,J)-DX1*P1(2,J)
          A1LM(J)=-DY2
          A2LM(J)=DX1
          A3LM(J)=0.
          XLIMS1(1,J)=MIN(P1(1,J),P2(1,J))
          XLIMS2(1,J)=MAX(P1(1,J),P2(1,J))
          YLIMS1(1,J)=MIN(P1(2,J),P2(2,J))
          YLIMS2(1,J)=MAX(P1(2,J),P2(2,J))
          IF (P1(1,J).EQ.P2(1,J)) THEN
            XLIMS1(1,J)=XLIMS1(1,J)-0.1
            XLIMS2(1,J)=XLIMS2(1,J)+0.1
          ENDIF
          IF (P1(2,J).EQ.P2(2,J)) THEN
            YLIMS1(1,J)=YLIMS1(1,J)-0.1
            YLIMS2(1,J)=YLIMS2(1,J)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND Z-INDEPENDENCE OF SURFACE EQUATION
          P3(3,J)=1.D55
          DL=SQRT(DX1**2+DY2**2)
          IF (NLTRZ) THEN
            DZ=ZLIMS2(1,J)-ZLIMS1(1,J)
            IF (DZ.GT.1.D20) DZ=ZDF
          ELSEIF (NLTRA) THEN
            XS=(P1(1,J)+P2(1,J))*0.5+RMTOR
            IF (ILTOR(J).GT.0) THEN
              DZ=ZLIMS2(1,J)-ZLIMS1(1,J)
C             IF (DZ.GT.1.D20) ????
            ELSEIF (ILTOR(J).EQ.0) THEN
              DZ=ZLIMS2(1,J)-ZLIMS1(1,J)
              IF (DZ.GT.1.D20) DZ=XS*TANAL/ALPHA*PI2A
            ENDIF
          ENDIF
          SAREA(J)=DL*DZ
          GOTO 90
        ELSEIF (RLB(J).GE.3.AND.RLB(J).LT.4.) THEN
C  1 TRIANGLE
          P4(1,J)=P3(1,J)
          P4(2,J)=P3(2,J)
          P4(3,J)=P3(3,J)
          P5(1,J)=P4(1,J)
          P5(2,J)=P4(2,J)
          P5(3,J)=P4(3,J)
        ELSEIF (RLB(J).GE.4..AND.RLB(J).LT.5) THEN
C  1 QUADRANGLE = 2 TRIANGLES
          P5(1,J)=P4(1,J)
          P5(2,J)=P4(2,J)
          P5(3,J)=P4(3,J)
        ENDIF
C
        TEST4=0.
        TEST5=0.
C
C  SET PLANE SURFACE FROM CO-ORDINATES OF THE FIRST 3 VERTICES
C  RLB.GE.3.0 AT THIS POINT
C
        DO 82 I=1,3
          XB(I)=P1(I,J)-P3(I,J)
          XC(I)=P2(I,J)-P3(I,J)
          XD(I)=P2(I,J)-P4(I,J)
          XE(I)=P3(I,J)-P4(I,J)
          XF(I)=P3(I,J)-P5(I,J)
          XG(I)=P4(I,J)-P5(I,J)
82      CONTINUE
C
        A1LM(J)=XB(2)*XC(3)-XB(3)*XC(2)
        A2LM(J)=XB(3)*XC(1)-XB(1)*XC(3)
        A3LM(J)=XB(1)*XC(2)-XB(2)*XC(1)
        A0LM(J)=-(A1LM(J)*P1(1,J)+A2LM(J)*P1(2,J)+A3LM(J)*P1(3,J))
C
C   SURFACE AREA: SAREA
        B=SQRT(XB(1)**2+XB(2)**2+XB(3)**2+EPS60)
        C=SQRT(XC(1)**2+XC(2)**2+XC(3)**2+EPS60)
        D=SQRT(XD(1)**2+XD(2)**2+XD(3)**2+EPS60)
        E=SQRT(XE(1)**2+XE(2)**2+XE(3)**2+EPS60)
        F=SQRT(XF(1)**2+XF(2)**2+XF(3)**2+EPS60)
        G=SQRT(XG(1)**2+XG(2)**2+XG(3)**2+EPS60)
        CG1=(XB(1)*XC(1)+XB(2)*XC(2)+XB(3)*XC(3))/B/C
        CG2=(XD(1)*XE(1)+XD(2)*XE(2)+XD(3)*XE(3))/D/E
        CG3=(XF(1)*XG(1)+XF(2)*XG(2)+XF(3)*XG(3))/F/G
        AR1=0.5*B*C*SQRT(1.-CG1*CG1+EPS60)
        AR2=0.5*D*E*SQRT(1.-CG2*CG2+EPS60)
        AR3=0.5*F*G*SQRT(1.-CG3*CG3+EPS60)
        SAREA(J)=AR1+AR2+AR3
C
C   TEST, IF 4TH  AND 5TH POINT ON SURFACE
        TEST4=A0LM(J)+A1LM(J)*P4(1,J)+A2LM(J)*P4(2,J)+A3LM(J)*P4(3,J)
        TEST5=A0LM(J)+A1LM(J)*P5(1,J)+A2LM(J)*P5(2,J)+A3LM(J)*P5(3,J)
C
C  NEW CO-ORDINATE SYSTEM IN P3,P1,P2 ; P4,P2,P3 ; P5,P3,P4
C    XS = P3 + XLS1*(P1-P3) + XMS1*(P2-P3)
C    XS = P4 + XLS2*(P2-P4) + XMS2*(P3-P4)
C    XS = P5 + XLS3*(P3-P5) + XMS3*(P4-P5)
C  PREPARE ARRAYS FOR COMPUTATION OF XLS1,...XMS3
        V1V23=XB(1)*XC(1)+XB(2)*XC(2)+XB(3)*XC(3)
        V113=XB(1)*XB(1)+XB(2)*XB(2)+XB(3)*XB(3)
        V223=XC(1)*XC(1)+XC(2)*XC(2)+XC(3)*XC(3)
C
        IF (ABS(V113).LE.EPS12) GOTO 98
        HELP=V1V23*V1V23/V113-V223
        IF (ABS(HELP).LE.EPS12) GOTO 98
        XK3=V1V23/V113
        XN3=1./HELP
        PS13(1,J)=(XB(1)*XK3-XC(1))*XN3
        PS13(2,J)=(XB(2)*XK3-XC(2))*XN3
        PS13(3,J)=(XB(3)*XK3-XC(3))*XN3
        P1A(J)=-P3(1,J)*PS13(1,J)-P3(2,J)*PS13(2,J)-P3(3,J)*PS13(3,J)
C
        IF (ABS(V223).LE.EPS12) GOTO 98
        HELP=V1V23*V1V23/V223-V113
        IF (ABS(HELP).LE.EPS12) GOTO 98
        XK3=V1V23/V223
        XN3=1./HELP
        PS23(1,J)=(XC(1)*XK3-XB(1))*XN3
        PS23(2,J)=(XC(2)*XK3-XB(2))*XN3
        PS23(3,J)=(XC(3)*XK3-XB(3))*XN3
        P2A(J)=-P3(1,J)*PS23(1,J)-P3(2,J)*PS23(2,J)-P3(3,J)*PS23(3,J)
C
        IF (RLB(J).GE.4) THEN
          V2V34=XD(1)*XE(1)+XD(2)*XE(2)+XD(3)*XE(3)
          V224=XD(1)*XD(1)+XD(2)*XD(2)+XD(3)*XD(3)
          V334=XE(1)*XE(1)+XE(2)*XE(2)+XE(3)*XE(3)
C
          IF (ABS(V224).LE.EPS12) GOTO 98
          HELP=V2V34*V2V34/V224-V334
          IF (ABS(HELP).LE.EPS12) GOTO 98
          XK4=V2V34/V224
          XN4=1./HELP
          PS24(1,J)=(XD(1)*XK4-XE(1))*XN4
          PS24(2,J)=(XD(2)*XK4-XE(2))*XN4
          PS24(3,J)=(XD(3)*XK4-XE(3))*XN4
          P1B(J)=-P4(1,J)*PS24(1,J)-P4(2,J)*PS24(2,J)-P4(3,J)*PS24(3,J)
C
          IF (ABS(V334).LE.EPS12) GOTO 98
          HELP=V2V34*V2V34/V334-V224
          IF (ABS(HELP).LE.EPS12) GOTO 98
          XK4=V2V34/V334
          XN4=1./HELP
          PS34(1,J)=(XE(1)*XK4-XD(1))*XN4
          PS34(2,J)=(XE(2)*XK4-XD(2))*XN4
          PS34(3,J)=(XE(3)*XK4-XD(3))*XN4
          P2B(J)=-P4(1,J)*PS34(1,J)-P4(2,J)*PS34(2,J)-P4(3,J)*PS34(3,J)
C
          IF (RLB(J).GE.5) THEN
            V3V45=XF(1)*XG(1)+XF(2)*XG(2)+XF(3)*XG(3)
            V335=XF(1)*XF(1)+XF(2)*XF(2)+XF(3)*XF(3)
            V445=XG(1)*XG(1)+XG(2)*XG(2)+XG(3)*XG(3)
C
            IF (ABS(V335).LE.EPS12) GOTO 98
            HELP=V3V45*V3V45/V335-V445
            IF (ABS(HELP).LE.EPS12) GOTO 98
            XK5=V3V45/V335
            XN5=1./HELP
            PS35(1,J)=(XF(1)*XK5-XG(1))*XN5
            PS35(2,J)=(XF(2)*XK5-XG(2))*XN5
            PS35(3,J)=(XF(3)*XK5-XG(3))*XN5
            P1C(J)=-P5(1,J)*PS35(1,J)-P5(2,J)*PS35(2,J)-
     -              P5(3,J)*PS35(3,J)
C
            IF (ABS(V445).LE.EPS12) GOTO 98
            HELP=V3V45*V3V45/V445-V335
            IF (ABS(HELP).LE.EPS12) GOTO 98
            XK5=V3V45/V445
            XN5=1./HELP
            PS45(1,J)=(XG(1)*XK5-XF(1))*XN5
            PS45(2,J)=(XG(2)*XK5-XF(2))*XN5
            PS45(3,J)=(XG(3)*XK5-XF(3))*XN5
            P2C(J)=-P5(1,J)*PS45(1,J)-P5(2,J)*PS45(2,J)-
     -              P5(3,J)*PS45(3,J)
          ENDIF
        ENDIF
C
        IF (ABS(TEST4).GT.EPS10) THEN
          WRITE (6,*) 'WARNING FROM TIMEA0, TEST FOR 4.TH POINT:'
          WRITE (6,*) 'SF. NUMBER, TEST= ',J,TEST4
        ELSEIF (ABS(TEST5).GT.EPS10) THEN
          WRITE (6,*) 'WARNING FROM TIMEA0, TEST FOR 5.TH POINT:'
          WRITE (6,*) 'SF. NUMBER, TEST= ',J,TEST5
        ENDIF
C
C  SET SOME MORE ASSISTENT SURFACE DATA AND CHECK CONSISTENCY
C  ST. NO. 90 --- 99
C
90      CONTINUE
C
        RLBNOT(J)=RLB(J).EQ.1.5.OR.RLB(J).EQ.2.5.OR.RLB(J).EQ.3.5.OR.
     .            RLB(J).EQ.4.5.OR.RLB(J).EQ.5.5
C
        IF (ILSWCH(J).NE.0.AND.ILIIN(J).EQ.0) THEN
          WRITE (6,*) 'EXIT FROM TIMEA0: SURFACE NO J IS OPERATING'
          WRITE (6,*) 'A SWITCH BUT IS NOT SEEN BY HISTORY'
          WRITE (6,*) 'J= ',J
          CALL EXIT
        ENDIF
C
        IF (ILSWCH(J).NE.0.AND.ILIIN(J).GT.0) THEN
          DO ISPZ=1,NSPTOT
            IF (TRANSP(ISPZ,1,J).NE.0.D0.OR.TRANSP(ISPZ,2,J).NE.0.D0)
     .      GOTO 95
          ENDDO
          GOTO 96
95        WRITE (6,*) 'EXIT FROM TIMEA0: SURFACE NO J IS OPERATING'
          WRITE (6,*) 'A SWITCH BUT IS SOMETIMES TRANSPARENT AND '
          WRITE (6,*) 'SOMETIMES REFLECTING (SEMI-TRANSPARENCY OPTION)'
          WRITE (6,*) 'POSSIBLE FIXES: MANUAL, CHAPTER 2, SECTION 6 '
          WRITE (6,*) 'J= ',J
          CALL EXIT
        ENDIF
C
96      GOTO 94
C
98      WRITE (6,*) 'WARNING FROM TIMEA0, SURFACE NO J IS TURNED OF'
        WRITE (6,*) 'BECAUSE ILL DEFINED '
        WRITE (6,*) 'J= ',J
        IGJUM0(J)=1
        IF (NLIMPB >= NLIMPS) THEN
          DO IAB=0,NLIMI
            IGJUM1(IAB,J)=1
          END DO
        ELSE
          DO IAB=0,NLIMI
            CALL BITSET (IGJUM1,0,NLIMPS,IAB,J,1,NBITS)
          END DO
        END IF
        GOTO 97
C
94      CONTINUE
C
        JUMLIM(J)=0
        ALM(J)=2.*A4LM(J)
        BLM(J)=2.*A5LM(J)
        CLM(J)=2.*A6LM(J)
        IF (A4LM(J).EQ.0..AND.A5LM(J).EQ.0..AND.A6LM(J).EQ.0..AND.
     .      A7LM(J).EQ.0..AND.A8LM(J).EQ.0..AND.A9LM(J).EQ.0.) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM1(J,J)=1
          ELSE
            CALL BITSET (IGJUM1,0,NLIMPS,J,J,1,NBITS)
          END IF
          AT=MAX(ABS(A1LM(J)),ABS(A2LM(J)),ABS(A3LM(J)))
          IF (ABS(A1LM(J)).EQ.AT) JUMLIM(J)=1
          IF (ABS(A2LM(J)).EQ.AT) JUMLIM(J)=2
          IF (ABS(A3LM(J)).EQ.AT) JUMLIM(J)=3
          XNORM=SQRT(A1LM(J)*A1LM(J)+A2LM(J)*A2LM(J)+A3LM(J)*A3LM(J))
          A0LM(J)=A0LM(J)/XNORM
          A1LM(J)=A1LM(J)/XNORM
          A2LM(J)=A2LM(J)/XNORM
          A3LM(J)=A3LM(J)/XNORM
          JUM=JUMLIM(J)
          GOTO (91,92,93),JUM
91          ALM(J)=-A0LM(J)/A1LM(J)
            BLM(J)=-A2LM(J)/A1LM(J)
            CLM(J)=-A3LM(J)/A1LM(J)
          GOTO 97
92          ALM(J)=-A0LM(J)/A2LM(J)
            BLM(J)=-A1LM(J)/A2LM(J)
            CLM(J)=-A3LM(J)/A2LM(J)
          GOTO 97
93          ALM(J)=-A0LM(J)/A3LM(J)
            BLM(J)=-A1LM(J)/A3LM(J)
            CLM(J)=-A2LM(J)/A3LM(J)
        ENDIF
97    CONTINUE
C
      CALL LEER(2)
99    CONTINUE
      RETURN
C
      ENTRY TIMEA1(MSURF,NCELL,NLI,NLE,NTCELL,IPERID,XX,YY,ZZ,TMT,
     .             VXX,VYY,VZZ,VV,
     .             MASURF,XR,YR,ZR,SG,TL,NLTRC)
      LM1=NLI
      LM2=NLE
C  PARTICLE ON STANDARD SURFACE?
      IF (MSURF.GT.NLIM) MSURF=0
C  FIND LOCAL COORDINATE SYSTEM IN CASE OF NLTRA
      IF (NLTRA.AND.NLTOR) THEN
        NNTCL=NTCELL
      ELSEIF (NLTRA.AND..NOT.NLTOR) THEN
        NNTCL=IPERID
      ENDIF
C  SAVE INITIAL CO-ORDINATES
      XS=XX
      YS=YY
      ZS=ZZ
      TS=TMT
      VXS=VXX
      VYS=VYY
      VZS=VZZ
      VVS=VV
      NNTCLS=NNTCL
C  WORKING CO-ORDINATES
      X=XX
      Y=YY
      Z=ZZ
      T=TMT
      VX=VXX
      VY=VYY
      VZ=VZZ
      V=VV
      NN=NNTCL
c slmod begin (sl)
c...  Suppress motion in the z-direction:
      IF ( OPTZMOTION.EQ.1.OR.
     .    (OPTZMOTION.EQ.2.AND.NACELL1.EQ.0).OR.
     .    (OPTZMOTION.EQ.3.AND.(NACELL1.EQ.0.OR.
     .     (XX.LT.62.5.AND.YY.GT.-61.0))).OR.
     .    (OPTZMOTION.EQ.4.AND.NACELL1.GT.0.AND.
     .     (XX.GT.62.5.OR.YY.LT.-61.0)).OR.
     .    (OPTZMOTION.EQ.5.AND.NACELL1.GT.0)) THEN
        VZS=0.0D0
        VZ=0.0D0
      ENDIF
c slmod end
C
      NLLLI=MSURF
c      IF (NLTRC) THEN
c        CALL LEER(1)
c        WRITE (6,*) 'TIMEA ,X,Y,Z,T ',X,Y,Z,T
c        WRITE (6,*) '      VX,VY,VZ,V ',VX,VY,VZ,V
c        IF (NLTRA) WRITE (6,*) 'MSURF,NNTCL ',MSURF,NNTCL
c        IF (.NOT.NLTRA) WRITE (6,*) 'MSURF ',MSURF
c      ENDIF
C
      TADD=0.
C
1000  TMIN=1.D30
      TL=1.D30
      MASURF=0
c slmod begin 
c...  Not sure how useful this is to other EIRENE users?  LM1 is reset so 
c     collisons are only calculated for a subset of additional surfaces,
c     which can save a lot of time if NLIMI is large. Code indent
c     is not used.  I should replace SBGKI, EBGKI, etc. with something
c     more generic:
      DO IREG=1,IGJUM4(NCELL,0)
        IF     (IGJUM4(NCELL,IREG).EQ. 0) THEN
        ELSEIF (IGJUM4(NCELL,IREG).EQ.-1) THEN
          LM1=1			
          LM2=SBGKI-1		
        ELSEIF (IGJUM4(NCELL,IREG).EQ.-2) THEN
          LM1=EBGKI+1
          LM2=HADDI
        ELSEIF (IGJUM4(NCELL,IREG).EQ.-3) THEN
          LM1=HADDI+1
          LM2=HSTDI
        ELSEIF (IGJUM4(NCELL,IREG).EQ.-4) THEN
          LM1=HSTDI+1
          LM2=NLIMI
        ELSE
          LM1=IGJUM4(NCELL,IREG)
          LM2=LM1
        ENDIF
c slmod end
C
C  LOOP OVER SURFACE NUMBER, DO 100
C
      DO 100 J=LM1,LM2
        IF (IGJUM0(J).NE.0) GOTO 100
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM1(NLLLI,J) .NE. 0) GOTO 100
        ELSE
          IF (BITGET(IGJUM1,0,NLIMPS,NLLLI,J,NBITS)) GOTO 100
        END IF
        IF (NCELL.LE.NOPTIM) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IF (IGJUM3(NCELL,J).NE.0) GOTO 100
          ELSE
            IF (BITGET(IGJUM3,0,NOPTIM,NCELL,J,NBITS)) GOTO 100
          END IF
        ENDIF
C
C  FOR NLTRA OPTION ONLY:
C  X,Z,VX AND VZ ARE GIVEN IN TOROIDAL CELL NN,
C  TRANSFORM CO-ORDINATES FOR THIS TRACK FROM LOCAL SYSTEM NN
C  TO THE LOCAL SYSTEM ILTOR(J), IN WHICH SURFACE J IS GIVEN
C  IF (ILTOR(J).LE.0) THIS SURFACE HAS TOROIDAL SYMMETRY
C
        IF (NLTRA) THEN
c slmod begin (sl)
c...      Toroidal resolution outside the standard grid.  Replace
c         this option check, with something custom.  The surfaces
c         for a particular toroidal segment *never* need to be rotated
c         here, because they already are in DIVIMP:
          IF (.NOT.NLTOR) THEN
            IF (ILTOR(J).EQ.0.OR.ILTOR(J).EQ.NTRSEG) THEN
              NN=NNTCLS
              X=XS
              Z=ZS
              VX=VXS
              VZ=VZS
            ELSE            
              GOTO 100
            ENDIF
          ELSEIF (ILTOR(J).GT.0.AND.ILTOR(J).NE.NN) THEN
c
c          IF (ILTOR(J).GT.0.AND.ILTOR(J).NE.NN) THEN
c slmod end
            CALL FZRTOR (X,Z,NN,XXR,PPP,NTNEW,.FALSE.,0)
            CALL FZRTRI (X,Z,ILTOR(J),XXR,PPP,NTNEW)
            ROT=2.*(NN-ILTOR(J))*ALPHA
            VSAVE=VX
            VX=COS(ROT)*VSAVE-SIN(ROT)*VZ
            VZ=SIN(ROT)*VSAVE+COS(ROT)*VZ
            NN=ILTOR(J)
C           WRITE (6,*) 'J,ILTOR(J),X,Z,VX,VZ ',J,ILTOR(J),X,Z,VX,VZ
C  X,Z,VX AND VZ ARE NOW GIVEN IN CELL NN=ILTOR(J). SO ARE THE COEFFICIENTS
C  OF SURFACE NO. J. FIND INTERSECTION IN THIS LOCAL SYSTEM
          ELSEIF (ILTOR(J).EQ.0) THEN
C  TOROIDALLY SYMMETRIC SURFACE, SURFACE COEFFICIENTS ARE THE SAME
C  IN EACH TOROIDAL CELL, THUS ESPECIALLY IN CELL NN
            NN=NNTCLS
            X=XS
            Z=ZS
            VX=VXS
            VZ=VZS
C           WRITE (6,*) 'J,ILTOR(J),X,Z,VX,VZ ',J,ILTOR(J),X,Z,VX,VZ
          ENDIF
        ENDIF
C
C  FIND INTERSECTION TIME TMX WITH BOUNDARY NO. J
C
        A1=0.
        GOTO (60,63,66),JUMLIM(J)
C  A1*TMX*TMX+A2*TMX+A3=0
        A1=(A4LM(J)*VX+A7LM(J)*VY+A8LM(J)*VZ)*VX+
     .     (A5LM(J)*VY+A9LM(J)*VZ)*VY+A6LM(J)*VZ*VZ
        A2=(A1LM(J)+ALM(J)*X)*VX+(A2LM(J)+BLM(J)*Y)*VY+
     .     (A3LM(J)+CLM(J)*Z)*VZ+
     .      A7LM(J)*(VX*Y+VY*X)+A8LM(J)*(VX*Z+VZ*X)+A9LM(J)*(VY*Z+VZ*Y)
C  A3 =  0. ?
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM2(NLLLI,J).NE.0) GOTO 40
        ELSE
          IF (BITGET(IGJUM2,0,NLIMPS,NLLLI,J,NBITS)) GOTO 40
        END IF
C  NO
        A3=A0LM(J)+(A1LM(J)+A4LM(J)*X+A7LM(J)*Y+A8LM(J)*Z)*X
     .            +(A2LM(J)+A5LM(J)*Y+A9LM(J)*Z)*Y
     .            +(A3LM(J)+A6LM(J)*Z)*Z
        IF (A1.EQ.0.) GOTO 50
        F=-A2/(A1+A1)
        G=F*F-A3/A1
        IF (G.LT.0.) GOTO 100
        G=SQRT(G)
        TMA=F+G
        TMI=F-G
c        IF (NLTRC) WRITE (6,*) 'TIMEA, J,TMI,TMA ',J,TMI,TMA
        IF (TMA.LE.EPS12.OR.TMI.GT.TMIN) GOTO 100
        IF (RLB(J).LT.0.) GOTO 21
        IF (RLB(J).GT.0.) GOTO 31
C
C  RLB(J) .EQ. 0.  STATEMENT 11---20
C
11      IF (TMI.LE.EPS12) GOTO 12
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        TMX=TMI
        GOTO 70
C
12      IF (TMA.GT.TMIN) GOTO 100
        WR=A2+A1*(TMA+TMA)
16      XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
18      TMX=TMA
        GOTO 70
C
C  CHECK BOUNDARY INEQUALITIES OF SURFACE
C  LGJ=.TRUE.: INTERSECTION POINT IS STILL VALID
C  LGJ=.FALSE.: INTERSECTION POINT IS OUTSIDE THE SPECIFIED AREA
C
C  RLB(J) .LT. 0.  STATEMENT 21---30
C
21      IF (TMI.LE.EPS12) GOTO 22
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        TUP=TMI
        ICOUNT=1
        GOTO 28
C
22      IF (TMA.GT.TMIN) GOTO 100
        WR=A2+A1*(TMA+TMA)
26      XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
        TUP=TMA
        ICOUNT=2
C
28      LGJ=.TRUE.
        IF (ILIN(J).GT.0) THEN
          I=0
27        I=1+I
          TST=ALIMS(I,J)+XLIMS(I,J)*XN+YLIMS(I,J)*YN+ZLIMS(I,J)*ZN
          LGJ=TST.LE.0.
          IF (LGJ.AND.I.LT.ILIN(J)) GOTO 27
        ENDIF
        IF (LGJ.AND.ISCN(J).GT.0) THEN
          I=0
29        I=1+I
          TST=ALIMS0(I,J)+
     .        XN*(XLIMS1(I,J)+XN*XLIMS2(I,J)+YN*XLIMS3(I,J))+
     .        YN*(YLIMS1(I,J)+YN*YLIMS2(I,J)+ZN*ZLIMS3(I,J))+
     .        ZN*(ZLIMS1(I,J)+ZN*ZLIMS2(I,J)+XN*YLIMS3(I,J))
          LGJ=TST.LE.0.
          IF (LGJ.AND.I.LT.ISCN(J)) GOTO 29
        ENDIF
        IF (LGJ) THEN
          TMX=TUP
          GOTO 70
        ELSEIF (ICOUNT.EQ.1) THEN
          GOTO 22
        ENDIF
        GOTO 100
C
C   RLB(J) .GT. 0.  STATEMENT 31---40
C
31      IF (TMI.LE.EPS12) GOTO 32
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        LGJ=XN.LE.XLIMS2(1,J).AND.XN.GE.XLIMS1(1,J).AND.
     .      YN.LE.YLIMS2(1,J).AND.YN.GE.YLIMS1(1,J).AND.
     .      ZN.LE.ZLIMS2(1,J).AND.ZN.GE.ZLIMS1(1,J)
        IF (RLBNOT(J)) LGJ=.NOT.LGJ
        TMX=TMI
        IF (LGJ) GOTO 70
C
32      IF (TMA.GT.TMIN) GOTO 100
        WR=A2+A1*(TMA+TMA)
36      XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
38      IF (RLB(J).LT.2.) THEN
          LGJ=XN.LE.XLIMS2(1,J).AND.XN.GE.XLIMS1(1,J).AND.
     .        YN.LE.YLIMS2(1,J).AND.YN.GE.YLIMS1(1,J).AND.
     .        ZN.LE.ZLIMS2(1,J).AND.ZN.GE.ZLIMS1(1,J)
        ELSE
          XMS1=XN*PS13(1,J)+YN*PS13(2,J)+ZN*PS13(3,J)+P1A(J)
          XLS1=XN*PS23(1,J)+YN*PS23(2,J)+ZN*PS23(3,J)+P2A(J)
          LGJ=XMS1.GE.0..AND.XLS1.GE.0..AND.XMS1+XLS1.LE.1.
          IF (RLB(J).GE.4.AND..NOT.LGJ) THEN
            XMS2=XN*PS24(1,J)+YN*PS24(2,J)+ZN*PS24(3,J)+P1B(J)
            XLS2=XN*PS34(1,J)+YN*PS34(2,J)+ZN*PS34(3,J)+P2B(J)
            LGJ=XMS2.GE.0..AND.XLS2.GE.0..AND.XMS2+XLS2.LE.1.
            IF (RLB(J).GE.5.AND..NOT.LGJ) THEN
              XMS3=XN*PS35(1,J)+YN*PS35(2,J)+ZN*PS35(3,J)+P1C(J)
              XLS3=XN*PS45(1,J)+YN*PS45(2,J)+ZN*PS45(3,J)+P2C(J)
              LGJ=XMS3.GE.0..AND.XLS3.GE.0..AND.XMS3+XLS3.LE.1.
            ENDIF
          ENDIF
        ENDIF
        IF (RLBNOT(J)) LGJ=.NOT.LGJ
        TMX=TMA
        IF (LGJ) GOTO 70
        GOTO 100
C
C   A1*TMX+A2=0
C
40      TMA=-A2/A1
c        IF (NLTRC) WRITE (6,*) 'TIMEA AT 40, J,TMA,A1,A2 ',J,TMA,A1,A2
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        WR=-A2
        IF (RLB(J)) 26,16,36
C
C   A2*TMX+A3=0
C
50      IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
c        IF (NLTRC) WRITE (6,*) 'TIMEA AT 50, J,TMA,A2,A3 ',J,TMA,A2,A3
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        WR=A2
        IF (RLB(J)) 26,16,36
C
C  A1LM(J).NE.0
C
60      CONTINUE
        A2=A1LM(J)*VX+A2LM(J)*VY+A3LM(J)*VZ
        A3=A0LM(J)+A1LM(J)*X+A2LM(J)*Y+A3LM(J)*Z
        IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        WR=A2
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
        XN=ALM(J)+BLM(J)*YN+CLM(J)*ZN
        TUP=TMA
        ICOUNT=2
        IF (RLB(J)) 28,18,38
C
C  A2LM(J).NE.0
C
63      CONTINUE
        A2=A1LM(J)*VX+A2LM(J)*VY+A3LM(J)*VZ
        A3=A0LM(J)+A1LM(J)*X+A2LM(J)*Y+A3LM(J)*Z
        IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        WR=A2
        XN=X+TMA*VX
        ZN=Z+TMA*VZ
        YN=ALM(J)+BLM(J)*XN+CLM(J)*ZN
        TUP=TMA
        ICOUNT=2
        IF (RLB(J)) 28,18,38
C
C  A3LM(J).NE.0
C
66      CONTINUE
        A2=A1LM(J)*VX+A2LM(J)*VY+A3LM(J)*VZ
        A3=A0LM(J)+A1LM(J)*X+A2LM(J)*Y+A3LM(J)*Z
        IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        WR=A2
        XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=ALM(J)+BLM(J)*XN+CLM(J)*YN
        TUP=TMA
        ICOUNT=2
        IF (RLB(J)) 28,18,38
C
C
C   DATA RETURNED TO CALLING PROGRAM
C   TENTATIVELY FOR SURFACE NO. J
C
70      TL=TMX+TADD
        TMIN=TMX
        XR=XN
        YR=YN
        ZR=ZN
        NNR=NN
        VXJ=VX
        VZJ=VZ
        MASURF=J
        SG=SIGN(1._DP,WR)
c        IF (NLTRC) THEN
c          WRITE (6,*) 'TIMEA, TL,XR,YR,ZR,NNR,MASURF,SG '
c          WRITE (6,*)         TL,XR,YR,ZR,NNR,MASURF,SG
c        ENDIF
C
C  LOOP OVER SURFACE-INDEX FINISHED
C
100   CONTINUE
c slmod begin
c...  The ENDDO for the IGJUM4 loop:
      ENDDO
c slmod end
C
C **********************************************************************
C
c slmod begin (debug)
      IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10) THEN
        IF (MASURF.EQ.0) THEN
          WRITE(6,'(4X,A)') 'TIMEA1: NO INTERSECTION FOUND'
        ELSEIF (ILIIN(MASURF).GT.0) THEN
          WRITE(6,'(4X,A)') 'TIMEA1: NON-TRANSPARENT SURFACE'
        ELSEIF (ILSIDE(MASURF)*SG.LT.0) THEN
          WRITE(6,'(4X,A)') 'TIMEA1: TRANSPARENT BUT WRONG SIDE'
        ENDIF
      ENDIF
c slmod end
C  IF NO INTERSECTION FOUND, RETURN
      IF (MASURF.EQ.0) RETURN
C  INTERSECTION AT SURFACE NO. MASURF
C  IF NOT TRANSPARENT, RETURN
      IF (ILIIN(MASURF).GT.0) RETURN
C  IF TRANSPARENT BUT WRONG SIDE, RETURN
      IF (ILSIDE(MASURF)*SG.LT.0) RETURN
C
c      IF (NLTRC) WRITE (6,*) 'NOT RETURNED FROM TIMEA, OTHER LOOP '
C  ILIIN=0, CONTINUE WITH ANOTHER LOOP IN SUBR. TIMEA
      IF (ILIIN(MASURF).EQ.0) THEN
        X=XR
        XS=X
        Y=YR
        YS=Y
        Z=ZR
        ZS=Z
C
        TADD=TL
        NLLLI=MASURF
        NN=NNR
        NNTCLS=NN
        VX=VXJ
        VXS=VXJ
        VZ=VZJ
        VZS=VZJ
        GOTO 1000
      ELSEIF (ILIIN(MASURF).LT.0) THEN
C  TRANSPARENT, BUT SWITCH AND/OR SURFACE TALLIES
c slmod begin (debug)
        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(4X,A)') 'TIMEA1: TRANSPARENT BUT SWITCH'
          WRITE(6,'(4X,A,I5,E12.5)')
     .      'TIMEA1: Output (MASURF,TL) ',
     .      masurf,tl
          WRITE(6,'(4X,A)') 'TIMEA1: '
        ENDIF
c slmod end
        RETURN
      ENDIF
C
      END
c === ROUTINE: timep
C
C
      SUBROUTINE TIMEP (ZRAD)
C
C   CALCULATE TIME SEGMENTS IN Y- OR POLOIDAL MESH
C
C   INPUT:
C
C       ZRAD   =
C       NRCELL = RADIAL CELL NUMBER, IN WHICH THIS TRACK OF LENGTH ZRAD
C                                    IS PERFORMED
C              = 0, IF TRACK OUTSIDE RADIAL MESH
C                   IF MRSURF .NE. 0 THEN
C                   REENTRY FOUND AT RADIAL SURFACE MRSURF.
C                   OTHERWISE: REENTRY AT NON DEFAULT
C                   POLOIDAL SURFACES IS SEARCHED FOR IN THIS CALL
C       NTCELL =
C       ITCELL =
C       NPCELL = POLOIDAL CELL INDEX OF POINT X00,Y00,Z00
C                AT WHICH THIS TRACK OF LENGTH ZRAD STARTS
C  WARNING: IF NLSRFY, NPCELL MAY BE WRONG
C       IPCELL =
C       X00,Y00,Z00: STARTING POINT FOR THIS TRACK
C
C       NCOUT,BLPD,NCOUNT
C
C   OUTPUT:
C
C       NCOU   = TOTAL NUMBER OF CELLS, WHICH THE CURRENT TRACK OF
C                LENGTH ZRAD HAS CROSSED
C       CLPD(I)= LENGTH OF THE PART NO. I (I=1,NCOU)
C       JUPC(I)= CELL NUMBER IN P-GRID FOR TRACK I
C       LUPC(I)= SURFACE NUMBER IN P-GRID AN THE END OF TRACK I
C       MUPC(I)= ORIENTATION OF TRACK I
C       NPCELL = LAST POLOIDAL CELL INDEX OF X00,Y00,Z00 POINT,
C                ON WHICH THIS TRACK OF LENGTH ZRAD ENDS
C       IPCELL = NUMBER OF FINAL POLOIDAL CELL NPCELL
C                IF TRACK ORIGINATED INSIDE STANDARD MESH
C              = NUMBER OF POLOIDAL CELL, AT WHICH REENTRY WAS FOUND
C                ON RADIAL OR TOROIDAL SURFACE
C                IF TRACK ORIGINATED OUTSIDE  STANDARD MESH
C       IRCELL = NUMBER OF RADIAL CELL, AT WHICH REENTRY AT POLOIDAL
C                SURFACE MPSURF WAS FOUND. OTHERWISE: UNMODIFIED.
C
C              = 0  IF TRACK COMPLETELY OUSIDE STANDARD MESH
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE COMPRT
      USE COMSOU
      USE CLGIN

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: ZRAD
c slmod begin
c...  Temp (DUM is temporarily declared in PARMMOD.f until the generalized
c     geometry code is integrated and tested):
      REAL(DP) :: HELP, GS, GC, F, XNEN, SUM, BB, DY, X0TEST,
c
c      REAL(DP) :: HELP, GS, GC, F, XNEN, DUM, SUM, BB, DY, X0TEST,
c slmod end
     .          Y0TEST, V1, V2, T1, ZRADS, ZRD, WIN1, X000, Z0T,
     .          TWIN1, XT, Y000, TSS, X0T, Y0T
      INTEGER :: ILLZ, IANP, I, IENP, JHELP, ISW, JJC, J1, NRCLLP,
     .           LHELP, J2, INCY, JN, ICOU, NYSAVE, NHELP,
     .           LEARCA, I1, MPTEST, IR, IRSAVE, NCOUPE, ISTS, IADD,
     .           ICOUT, NCPEN, IPOLGS, NCPAN, J, LEARC2, NJC, LEARC1,
     .           IST, ITEST, JSH, NN ,IN
      INTEGER :: NCOUNS(N2ND+N3RD)
      LOGICAL :: LCUTY(N2NDPLG), LCUTX(N1ST), lnincz
      SAVE
C
      IF (NLTRC) THEN
        CALL LEER(1)
        IF (NRCELL.GT.0) THEN
          WRITE (6,*) 'TIMEP FROM INSIDE, NPANU ', NPANU
          WRITE (6,*) 'ZRAD,NRCELL,NPCELL '
          WRITE (6,*)  ZRAD,NRCELL,NPCELL
        ELSE
          WRITE (6,*) 'TIMEP FROM OUTSIDE, NPANU ', NPANU
          WRITE (6,*) 'MRSURF,MTSURF,ZRAD '
          WRITE (6,*)  MRSURF,MTSURF,ZRAD
        ENDIF
      ENDIF
C
      ICOUT=1
      IADD=0
      ZRADS=ZRAD
      IPOLGS=IPOLGN
      ncpan=0
      ncpen=0
C
      IF (NLTOR) THEN
C       NCOUT=NCOUT
        ZRAD=BLPD(1)
        NTCELL=NCOUNT(1)
        IF (NCOUT.GT.1) IPOLGN=0
        IF (NLTRC) WRITE (6,*) 'WG. NLTOR: ZRAD,NTCELL ',ZRAD,NTCELL
      ELSE
        NCOUT=1
C       ZRAD=ZRAD
C       NTCELL=1
      ENDIF
C
10000 CONTINUE
C
CDR NOV.99: ADDED BECAUSE OF FOLION OPTION, WITH DEFAULT B-FIELD (IN Z DIRECTION)
C
      IF (ABS(VELZ).EQ.1.D0) THEN
        NCOUP=1
        ALPD(NCOUP)=ZRAD
        JUPC(NCOUP)=NPCELL
C       X00=X00+ZRAD*VELX
C       Y00=Y00+ZRAD*VELY
        IPCELL=NPCELL
        MPSURF=0
        GOTO 5000
      ENDIF
C
      IF (LEVGEO.EQ.3) THEN
c slmod begin (debug)
        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(5X,A,5I4,3X,1P,E12.5)')
     .      'TIMEP: Input  (MPSURF,NR,NP,IR,IPOLGN ZRAD) ',
     .      MPSURF,NRCELL,NPCELL,IRCELL,IPOLGN,ZRAD
c slmod end
C
        IF (NRCELL.GT.0) GOTO 10
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT POLOIDAL SURFACES
C
        NCOUP=0
        ZRD=ZRAD
c slmod begin (debug)
        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(5X,A)') 'TIMEP: CHECKING FOR RE-ENTRY'
c slmod end
        DO 3 ISTS=1,NSTSI
          MPTEST=INUMP(ISTS,2)
          IF (MPTEST.NE.0) THEN
C  TEST POLOIDAL SURFACE NO. MPTEST FOR REENTRY
            DO 2 IR=1,NR1STM
              I1=IR+1
c slmod begin
              IF (GRIDOPT.EQ.1) THEN
                V1=(YVERT(IR,MPTEST,1)-Y00)*VELX-
     .             (XVERT(IR,MPTEST,1)-X00)*VELY
                V2=(YVERT(IR,MPTEST,2)-Y00)*VELX-
     .             (XVERT(IR,MPTEST,2)-X00)*VELY
c...            Temp:
                DUM=(YPOL(IR,MPTEST)-Y00)*VELX-
     .              (XPOL(IR,MPTEST)-X00)*VELY
                CALL CHECKNUM('ZAAH1',IR,MPTEST,V1,DUM)
                DUM=(YPOL(I1,MPTEST)-Y00)*VELX-
     .              (XPOL(I1,MPTEST)-X00)*VELY
                CALL CHECKNUM('ZAAH2',IR,MPTEST,V2,DUM)
              ELSE
                V1=(YPOL(IR,MPTEST)-Y00)*VELX-(XPOL(IR,MPTEST)-X00)*VELY
                V2=(YPOL(I1,MPTEST)-Y00)*VELX-(XPOL(I1,MPTEST)-X00)*VELY
              ENDIF
c
c              V1=(YPOL(IR,MPTEST)-Y00)*VELX-(XPOL(IR,MPTEST)-X00)*VELY
c              V2=(YPOL(I1,MPTEST)-Y00)*VELX-(XPOL(I1,MPTEST)-X00)*VELY
c slmod end
              LCUTX(IR)=V1*V2.LE.0.
2           CONTINUE
            DO 4 IR=1,NR1STM
              IF (LCUTX(IR)) THEN
c slmod begin
                IF (GRIDOPT.EQ.1) THEN
                  T1=((XVERT(IR,MPTEST,1)-X00)*VVTY(IR,MPTEST)-
     .                (YVERT(IR,MPTEST,1)-Y00)*VVTX(IR,MPTEST))
     .               /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)
c...              Test:
                  DUM=((XPOL(IR,MPTEST)-X00)*VVTY(IR,MPTEST)-
     .                 (YPOL(IR,MPTEST)-Y00)*VVTX(IR,MPTEST))
     .               /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)
                  CALL CHECKNUM('ZAAZ1',IR,MPTEST,T1,DUM)
                ELSE
                  T1=((XPOL(IR,MPTEST)-X00)*VVTY(IR,MPTEST)-
     .                (YPOL(IR,MPTEST)-Y00)*VVTX(IR,MPTEST))
     .               /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)
                ENDIF

c...            Ignore a *very* short collision time, since it is likley an inadvertant
c               registration of a collision with the surface that the particle is
c               currently on (is this possible?) -- turning this off for EIRENE02 (Mar 05, 2003):
                IF (.FALSE..AND.T1.LT.EPS12) THEN
                  WRITE(0,*) 'IGNORING POLOIDAL SURFACE COLLISION',NPANU
                  T1 = -1.0D0
                ENDIF
c
c                T1=((XPOL(IR,MPTEST)-X00)*VVTY(IR,MPTEST)-
c     .              (YPOL(IR,MPTEST)-Y00)*VVTX(IR,MPTEST))
c     .             /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)
c slmod end
c                IF (NLTRC) WRITE (6,*) 'IR,MPTEST,T1 ',IR,MPTEST,T1
                IF (T1.LT.0..OR.T1.GE.ZRD) GOTO 4
c                IF (NLTRC) WRITE (6,*) 'VALID INTERSECTION AT T1= ',T1
c slmod begin (debug)
                IF (printopt.GE.1.AND.printopt.LE.10)
     .            WRITE(6,'(5X,A,2I4,1P,E12.4)')
     .              'TIMEP: *** VALID INTERSECTION *** (IR,MPTEST T1) ',
     .              IR,MPTEST,T1
c slmod end
                NCOUP=1
                LUPC(NCOUP)=MPTEST
                IRSAVE=IR
                MUPC(NCOUP)=SIGN(1._DP,VELX*PPLNX(IR,MPTEST)+
     .                                 VELY*PPLNY(IR,MPTEST))
                JUPC(NCOUP)=1
                ZRD=T1
              ENDIF
4           CONTINUE
          ENDIF
3       CONTINUE
c slmod begin (sl)
        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(5X,A,3I4)')
     .      'TIMEP: NCOUP,LUPC,MPSURF=',NCOUP,LUPC(MAX(1,NCOUP)),MPSURF

c GERMANY: GETTING RID OF THE SECOND REQUIREMENT SO THAT PARTICLES
c          CAN BOUNCE OFF THE BACK OF THE OUTER TARGET SEQUENTIALLY.
c          I MAY BE ABLE TO GET AWAY WITH THIS BECAUSE OF THE
c          NEW CONSTRAINT ON T1 ABOVE (T1>1.0E-14 REQUIRED):

c... NEED TO ADD A QUALIFIER HERE, SO THAT THE ORIGINAL CONDITION IS PRESERVED
c    FOR STANDARD EIRENE RUNS -- USING GRIDOPT FOR NOW, BUT SHOULD BE SOMETHING
c    BETTER

c... IS THIS QUALIFIER PLACED CORRECTLY -- SHOULD BE .AND.?
        IF (NCOUP.GT.0.AND.
     .      (GRIDOPT.EQ.1.OR.LUPC(MAX(1,NCOUP)).NE.MPSURF)) THEN
c
c        IF (NCOUP.GT.0.AND.LUPC(MAX(1,NCOUP)).NE.MPSURF) THEN
c slmod end
C  REENTRY FOUND, REDUCE ZRAD TO T1
C  NCOUP=1 AT THIS POINT
          NCOUPE=1
          MPSURF=LUPC(NCOUPE)
          IPOLGN=LUPC(NCOUPE)
          NINCY=MUPC(NCOUPE)
          IRCELL=IRSAVE
          ZRAD=ZRD
          ISRFCL=0
          ALPD(NCOUPE)=ZRAD
          MRSURF=0
          MTSURF=0
          MASURF=0
          NINCX=0
          NINCZ=0
          IF (NLTRC) THEN
            WRITE (6,*) 'REENTRY FOUND, MPSURF,ZRAD = ',MPSURF,ZRAD
            WRITE (6,*) 'IRCELL ',IRCELL
          ENDIF
c slmod begin (debug)
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE (6,*) 'CHECK HERE FOR CORRECT POLOIDAL CELL NUMBER!'
c slmod end
          GOTO 31
        ELSE
C  NO REENTRY FOUND
          NCOUP=1
          JUPC(1)=1
          ALPD(1)=ZRAD
          IPCELL=IPOLGN
          MPSURF=0
          IF (NLTRC) THEN
            WRITE (6,*) 'NO REENTRY INTO POLOIDAL GRID FOUND '
          ENDIF
c slmod begin (debug)
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(5X,A,2F10.4)')
     .        'TIMEP: Updating position 5 (VELX,VELY) ',VELX,VELY
c slmod end
          X00=X00+ZRAD*VELX
          Y00=Y00+ZRAD*VELY
          GOTO 5000
        ENDIF
C
10      CONTINUE
C
C  PARTICLE INSIDE STANDARD MESH, RADIAL CELL NO. NRCELL
C
        IF (NCOUP.EQ.0) THEN
          WRITE (6,*) 'ERROR IN TIMEP: NCOUP=0'
c slmod begin (sl)
c...      Just trying to kill particle:
          WRITE(6,*) 'TIMEP: KILL PARTICLE.  NPANU = ',NPANU
          WRITE(0,*) 'TIMEP: KILL PARTICLE.  NPANU = ',NPANU
          ZRAD=-1.0D0
c          ZRAD=1.0D+30
c slmod end
          RETURN
        ENDIF
        IF (NLTRC) THEN
          WRITE (6,*) ' TIMEP IN NEIGHBOR PART '
          WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
          WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
          WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
          WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
        ENDIF
        NLSRFY=.FALSE.
        IF ((icout == 1) .and. 
     .    (SQRT((X0-X00)**2+(Y0-Y00)**2).GT.EPS10)) THEN
          DO J=1,NCOUP
            ALPD(J)=ALPD(J)-ZT
          ENDDO
        ENDIF
C
C  ACCOUNT FOR OTHER SURFACES INSIDE POLOIDAL MESH
C
        lnincz=.false.
        if (ityp==3) lnincz=(zrad < alpd(1))
        IF (ALPD(NCOUP).GT.ZRAD) THEN
          DO J=1,NCOUP
            IF (ALPD(J).GT.ZRAD) THEN
              ncpan=j+1
              ncpen=ncoup+1
              do jsh=ncoup,j+1,-1
                alpd(jsh+1)=alpd(jsh)-zrad
                jupc(jsh+1)=jupc(jsh)
                lupc(jsh+1)=lupc(jsh)
                mupc(jsh+1)=mupc(jsh)
              end do
              alpd(ncpan)=alpd(j)-zrad
              jupc(ncpan)=jupc(j)
              lupc(ncpan)=lupc(j)
              mupc(ncpan)=mupc(j)  
              ALPD(J)=ZRAD
              IPOLGN=JUPC(J)
              NCOUP=J
              LUPC(NCOUP)=0
              MUPC(NCOUP)=0
              GOTO 4711
            ENDIF
          ENDDO
4711      CONTINUE
        ENDIF
C
C   ADJUST ALPD AND ACCOUNT FOR "NONDEFAULT" POLOIDAL SURFACES
C
        TSS=0.
        IST=0
        DO J=1,NCOUP
          ALPD(J)=ALPD(J)-TSS
          TSS=TSS+ALPD(J)
C
          ITEST=INMP2I(NRCELL,LUPC(J),0)
          IN=ITEST+NLIM
!pb          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
          IF ((ityp==3).or.(ITEST.NE.0.AND.ILIIN(IN).NE.0)) THEN
C
C  TRACK ENDS ON ONE OF THE NON DEFAULT POLOIDAL SURFACES
C
            IF (NLTRC) THEN
              WRITE (6,*) ' TRACK TERMINATED'
              WRITE (6,*) ' ITEST,ILIIN ',ITEST,ILIIN(IN)
            ENDIF
            NCOUPE=J
            MPSURF=LUPC(NCOUPE)
            IPOLGN=LUPC(NCOUPE)
            NINCY=MUPC(NCOUPE)
            ZRAD=TSS
c slmod begin (debug)
            IF (printopt.GE.1.AND.printopt.LE.10)
     .        WRITE(6,'(5X,A,G12.5)')
     .          'TIMEP: *** TRACK ENDS ON NON-DEFAULT SURF *** (ZRAD)',
     .          zrad
c slmod end
            if ((ITEST.NE.0.AND.ILIIN(IN).NE.0).or. 
     .          (zrad < tl)) then
              ISRFCL=0
              MASURF=0
            end if
!pb
            if (lnincz) then         ! toroidal surface
              nincx=0
              nincy=0
              mrsurf=0
              mpsurf=0
!pb            elseif (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
!pb              NINCX=0
!pb              NINCZ=0
!pb              MRSURF=0
!pb              MTSURF=0
            elseif ((ityp==3).and.(nincx.ne.0)) then   ! radial surface
              nincy=0
              nincz=0
              mpsurf=0
              mtsurf=0
            else                     ! poloidal surface
              nincx=0
              nincz=0
              mrsurf=0
              mtsurf=0
            end if
!pb
            NCOUP=NCOUPE
            GOTO 311
          ENDIF
C
        ENDDO
C
C  LAST CELL, TRACK DOES NOT END ON A POLOIDAL SURFACE
C
C  INDEX IPOLGN OF LAST CELL NOT KNOWN ?
        IF (IPOLGN.EQ.0) THEN
c slmod begin (sl)
c...      Just trying to kill particle.  This problem is caused by inaccuracies
c         in the FindPoloidalCell routine.
          WRITE(0,*) 'TIMEP: ATTEMPTING TO KILL PARTICLE (AVOIDING'//
     .               ' LEARC1)'
          WRITE(6,*) 'TIMEP: ATTEMPTING TO KILL PARTICLE NPANU = ',NPANU
          ZRAD=-1.0D0
c          ZRAD=1.0D+30
          RETURN
c
c          X0T=X00+VELX*ZRAD
c          Y0T=Y00+VELY*ZRAD
c          NN=LEARC1(X0T,Y0T,Z0T,IPOLGN,
c     .              NRCELL,NRCELL,.FALSE.,.FALSE.,
c     .              NPANU,'TIMEP       ')
c slmod end
        ENDIF
C
        MPSURF=0
        NINCY=0
C
311     CONTINUE
        IF (NLTRC) THEN
          WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
          WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
          WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
          WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
        ENDIF
C
        X00=X00+ZRAD*VELX
        Y00=Y00+ZRAD*VELY
        NPCELL=JUPC(NCOUP)
        IPCELL=NPCELL
        GOTO 5000
C
C
30      CONTINUE
C
C  LAST CELL, TRACK DOES NOT END ON A POLOIDAL SURFACE
C
C  INDEX OF LAST CELL NOT KNOWN (E.G. DUE TO ADD. SURFACE) ?
CDR  ERROR: FALLS NLSRFY, GGFLS NPCELL FALSCH
        IF (IPOLGN.EQ.0) THEN
CDR       IF (NCOUP.EQ.0) THEN
CDR         IPOLGN=NPCELL
CDR       ELSE
c slmod begin (debug)
            IF (debugopt.NE.0) THEN
              WRITE(6,*)
              WRITE(0,'(5X,A)') 'TIMEP: ERROR: LEARC2 TIMEP CALL WILL'//
     .                          ' FAIL'
              WRITE(0,*) '   NPANU = ',npanu
            ENDIF
c slmod end
            X0T=X00+VELX*ZRAD
            Y0T=Y00+VELY*ZRAD
            IPOLGN=LEARC2(X0T,Y0T,NRCELL,NPANU,'TIMEP       ')
CDR       ENDIF
        ENDIF
C
        NCOUPE=NCOUP+1
        JUPC(NCOUPE)=IPOLGN
        LUPC(NCOUPE)=0
        MUPC(NCOUPE)=0
        ALPD(NCOUPE)=ZRAD-TSS
        MPSURF=0
        NINCY=0
C
31      CONTINUE
        NCOUP=NCOUPE
        IF (NLTRC) THEN
          WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
          WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
          WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
          WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
        ENDIF
C
        X00=X00+ZRAD*VELX
        Y00=Y00+ZRAD*VELY
        NPCELL=JUPC(NCOUP)
        IPCELL=NPCELL
        GOTO 5000
C
      ELSEIF (LEVGEO.EQ.2) THEN
        IF (NLCRC) THEN
C
C
C  DISTANCE TO NEXT X- OR RADIAL SURFACE KNOWN?
C
C       IF (TS.GE.1.D30) GOTO 1000
C
C  YES! ZRAD IS THE DISTANCE TRAVELLED IN X- OR RADIAL CELL NO. NRCELL
C
          NCOUP=1
          IF (NRCELL.LT.1.OR.NRCELL.GT.NR1STM) THEN
            X00=X00+VELX*ZRAD
            Y00=Y00+VELY*ZRAD
            WIN1=MOD(ATAN2(Y00,X00)+PI2A-PSURF(1),PI2A)+PSURF(1)
            NPCELL=WIN1/YDF*DBLE(NP2NDM)+1.
            GOTO 5000
          ENDIF
C
C  THE OLD POLOIDAL CELL INDEX IS: NPCELL
C  FIND THE NEW CELL INDEX : NJC
C
          X000=X00+VELX*ZRAD
          Y000=Y00+VELY*ZRAD
          WIN1=MOD(ATAN2(Y000,X000)+PI2A-PSURF(1),PI2A)+PSURF(1)
          NJC=WIN1/YDF*DBLE(NP2NDM)+1.
C
          TWIN1=0.
          NCOUP=0
          IF (NJC.EQ.NPCELL) GOTO 150
C
C   FIND ORIENTATION IN THETA-GRID
C
          XT=-Y00*VELX+X00*VELY
          NINCY=1
          IF (XT.LT.0.) NINCY=-1
C
C   CONTRIBUTION TO EACH THETA-CELL
C   NPCELL : STARTINDEX
C   JJC    : SURFACEINDEX
C   J1     : CELLINDEX
C   NJC    : ENDINDEX
          J1=NPCELL
100       JJC=J1
          IF (NINCY.EQ.1) JJC=JJC+1
C   TIMESTEP FROM X00,Y00 TO THETA-SURFACE, THETA=WIN
          GS=SINPH(JJC)
          GC=COSPH(JJC)
          XNEN=VELX*GS-VELY*GC
          F=(Y00*GC-X00*GS)/(XNEN+EPS60)
          NCOUP=NCOUP+1
          JUPC(NCOUP)=J1
          LUPC(NCOUP)=JJC
          ALPD(NCOUP)=F-TWIN1
          TWIN1=F
C
          J1=J1+NINCY
          IF (J1.EQ.0) J1=NP2NDM
          IF (J1.EQ.NP2ND) J1=1
!pb          IF (J1.NE.NJC) GOTO 100
          IF ((ityp.ne.3).and.(J1.NE.NJC)) GOTO 100
C
C   LAST THETA-CELL
C
150       NCOUP=NCOUP+1
          JUPC(NCOUP)=NJC
          ALPD(NCOUP)=ZRAD-TWIN1
          X00=X000
          Y00=Y000
          NPCELL=NJC
          MPSURF=0
          NINCY=0
          GOTO 5000
C
        ELSE
C
C  NEW PART: NOT NLCRC
C  PARTICLE INSIDE STANDARD MESH
C
          NCOUP=0
          NRCLLP=NRCELL+1
C
C   SEARCH FOR ALL POSSIBLE INTERSECTIONS WITHIN THE RADIAL CELL NRCELL
C
          DO 111 I=1,NP2ND
111         LCUTY(I)=.FALSE.
C
          DO 112 J=1,NP2ND
c slmod begin (sl)
c...        LEVGEO=2:
            IF (GRIDOPT.EQ.1) THEN
              WRITE(0,*) 'GRIDOPT: OPTION NOT SUPPORTED 003'
              STOP
            ENDIF
c slmod end
            V1=(YPOL(NRCELL,J)-Y00)*VELX-(XPOL(NRCELL,J)-X00)*VELY
            V2=(YPOL(NRCLLP,J)-Y00)*VELX-(XPOL(NRCLLP,J)-X00)*VELY
            LCUTY(J)=V1*V2.LE.0.
            IF (NLTRC) THEN
              IF (LCUTY(J)) WRITE (6,*) 'LCUTY=TRUE FOR ',J
            ENDIF
112       CONTINUE
          IF (NLSRFY) THEN
            LCUTY(MPSURF)=.FALSE.
            NLSRFY=.FALSE.
C  PSURF(1)=PSURF(NP2ND)
            IF (MPSURF.EQ.1) LCUTY(NP2ND)=.FALSE.
            IF (MPSURF.EQ.NP2ND) LCUTY(1)=.FALSE.
          ENDIF
C
          IANP=ILLZ(NP2ND,LCUTY,1)+1
          IENP=NP2ND-ILLZ(NP2ND,LCUTY,-1)
C
C   COMPUTE THE FLIGHT TIMES TO THE INTERSECTION POINTS
C
        DO 114 I=IANP,IENP
          IF (LCUTY(I)) THEN
c slmod begin (sl)
c...        LEVGEO=2:
            IF (GRIDOPT.EQ.1) THEN
              WRITE(0,*) 'GRIDOPT: OPTION NOT SUPPORTED 004'
              STOP
            ENDIF
c slmod end
            T1=((XPOL(NRCELL,I)-X00)*VVTY(NRCELL,I)-
     .          (YPOL(NRCELL,I)-Y00)*VVTX(NRCELL,I))
     .         /(VELX*VVTY(NRCELL,I)-VELY*VVTX(NRCELL,I)+EPS60)
            IF (NLTRC) WRITE (6,*) 'I,T1 ',I,T1
            IF (T1.LT.0..OR.T1.GE.ZRAD) GOTO 114
            IF (NLTRC) WRITE (6,*) 'VALID INTERSECTION AT T1= ',T1
            NCOUP=NCOUP+1
            LUPC(NCOUP)=I
            MUPC(NCOUP)=SIGN(1._DP,VELX*PPLNX(NRCELL,I)+
     .                             VELY*PPLNY(NRCELL,I))
            IF (MUPC(NCOUP).EQ.-1) THEN
              JUPC(NCOUP)=I
              ALPD(NCOUP)=T1
              IF (I.EQ.NP2ND) NCOUP=NCOUP-1
            ELSEIF (MUPC(NCOUP).EQ.1) THEN
              JUPC(NCOUP)=I-1
              ALPD(NCOUP)=T1
              IF (I.EQ.1) NCOUP=NCOUP-1
            ENDIF
          ENDIF
114     CONTINUE
C
C   REARRANGE THE FLIGHT TIMES IN ASCENDING ORDER
C
115     ISW=0
        DO 120 J=1,NCOUP-1
          IF (ALPD(J).GT.ALPD(J+1)) THEN
            ISW=ISW+1
            HELP=ALPD(J)
            ALPD(J)=ALPD(J+1)
            ALPD(J+1)=HELP
            JHELP=JUPC(J)
            JUPC(J)=JUPC(J+1)
            JUPC(J+1)=JHELP
            LHELP=LUPC(J)
            LUPC(J)=LUPC(J+1)
            LUPC(J+1)=LHELP
            NHELP=MUPC(J)
            MUPC(J)=MUPC(J+1)
            MUPC(J+1)=NHELP
          ENDIF
120     CONTINUE
        IF (ISW.GT.0.AND.NCOUP.GT.2) GOTO 115
        IF (NLTRC.AND.NCOUP.GT.0) THEN
          WRITE (6,*) ' NACH SORTIEREN '
          WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
          WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
          WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
          WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
        ENDIF
C
        DO 125 J=1,NCOUP-1
          IF (ABS(ALPD(J+1)-ALPD(J)).LE.EPS30) THEN
            IF (JUPC(J).LE.0.OR.JUPC(J).GE.NP2ND) THEN
              IF (NLTRC) THEN
                WRITE (6,*) ' VERTAUSCHE ALPD(',J,') UND (',J+1,')'
                WRITE (6,*) ' ALPD = ',ALPD(J),ALPD(J+1)
                WRITE (6,*) ' JUPC = ',JUPC(J),JUPC(J+1)
              ENDIF
              HELP=ALPD(J)
              ALPD(J)=ALPD(J+1)
              ALPD(J+1)=HELP
              JHELP=JUPC(J)
              JUPC(J)=JUPC(J+1)
              JUPC(J+1)=JHELP
              LHELP=LUPC(J)
              LUPC(J)=LUPC(J+1)
              LUPC(J+1)=LHELP
              NHELP=MUPC(J)
              MUPC(J)=MUPC(J+1)
              MUPC(J+1)=NHELP
            ENDIF
          ENDIF
125     CONTINUE
C
C   ADJUST ALPD AND ACCOUNT FOR "NONDEFAULT" POLOIDAL SURFACES
C
        TSS=0.
        IST=0
        DO 130 J=1,NCOUP
          ALPD(J)=ALPD(J)-TSS
          TSS=TSS+ALPD(J)
C
          ITEST=INMP2I(NRCELL,LUPC(J),0)
          IN=ITEST+NLIM
!pb          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
          IF ((ityp==3).or.(ITEST.NE.0.AND.ILIIN(IN).NE.0)) THEN
C
C  TRACK ENDS ON ONE OF THE NON DEFAULT POLOIDAL SURFACES
C
            IF (NLTRC) THEN
              WRITE (6,*) ' TRACK TERMINATED'
              WRITE (6,*) ' ITEST,ILIIN ',ITEST,ILIIN(IN)
            ENDIF
            NCOUPE=J
            MPSURF=LUPC(NCOUPE)
            IPOLGN=LUPC(NCOUPE)
            NINCY=MUPC(NCOUPE)
            ZRAD=TSS
            ISRFCL=0
            NINCX=0
            NINCZ=0
            MRSURF=0
            MTSURF=0
            MASURF=0
            GOTO 131
          ENDIF
C
130     CONTINUE
C
C  LAST CELL, TRACK DOES NOT END ON A POLOIDAL SURFACE
C
C  INDEX OF LAST CELL NOT KNOWN (E.G. DUE TO ADD. SURFACE) ?
CDR  ERROR: FALLS NLSRFY, GGFLS NPCELL FALSCH
        IF (IPOLGN.EQ.0) THEN
CDR       IF (NCOUP.EQ.0) THEN
CDR         IPOLGN=NPCELL
CDR       ELSE
            X0T=X00+VELX*ZRAD
            Y0T=Y00+VELY*ZRAD
            IPOLGN=LEARC2(X0T,Y0T,NRCELL,NPANU,'TIMEP       ')
CDR       ENDIF
        ENDIF
C
        NCOUPE=NCOUP+1
        JUPC(NCOUPE)=IPOLGN
        LUPC(NCOUPE)=0
        MUPC(NCOUPE)=0
        ALPD(NCOUPE)=ZRAD-TSS
        MPSURF=0
        NINCY=0
C
131       CONTINUE
          NCOUP=NCOUPE
          IF (NLTRC) THEN
            WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
            WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
            WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
            WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
          ENDIF
C
c slmod begin (debug)
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(5X,A,2F10.4)')
     .        'TIMEP: Updating position 7 (VELX,VELY) ',VELX,VELY
c slmod end
          X00=X00+ZRAD*VELX
          Y00=Y00+ZRAD*VELY
          NPCELL=JUPC(NCOUP)
          IPCELL=NPCELL
          GOTO 5000
C
        ENDIF
C
      ELSEIF (LEVGEO.EQ.1) THEN
C
C  IDENTICAL, UP TO NAMES, TO TOROIDAL GRID PART, NLTRZ BLOCK
C
        IF (NRCELL.GT.0) GOTO 2900
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT POLOIDAL SURFACES
C
        NCOUP=0
        ZRD=ZRAD
        BB=VELY+EPS60
        NYSAVE=1
        IF (VELY.LT.0.D0) NYSAVE=-1
        DO 2903 ISTS=1,NSTSI
          MPTEST=INUMP(ISTS,2)
          IF (MPTEST.NE.0) THEN
C  TEST POLOIDAL SURFACE NO. MPTEST FOR REENTRY
C  TIME FROM Y00 TO PSURF
            DY=PSURF(MPTEST)-Y00
            F=DY/BB
            IF (NLTRC) WRITE (6,*) 'MPTEST,F,DY ',MPTEST,F,DY
            IF (F.LE.ZRD.AND.F.GT.0.D0) THEN
              X0TEST=X00+VELX*F
              IF (X0TEST.GE.RSURF(1).AND.X0TEST.LE.RSURF(NR1ST)) THEN
                IRSAVE=LEARCA(X0TEST,RSURF,1,NR1ST,1,'TIMEP 1    ')
                NCOUP=1
                JUPC(NCOUP)=1
                LUPC(NCOUP)=MPTEST
                MUPC(NCOUP)=NYSAVE
                ZRD=F
              ENDIF
            ENDIF
          ENDIF
2903    CONTINUE
        IF (NCOUP.GT.0.AND.LUPC(MAX(1,NCOUP)).NE.MPSURF) THEN
C  REENTRY FOUND, REDUCE ZRAD TO F
C  NCOUP=1 AT THIS POINT
          NCOUP=1
          MPSURF=LUPC(NCOUP)
          NINCY=MUPC(NCOUP)
          IRCELL=IRSAVE
          ZRAD=ZRD
          ISRFCL=0
          ALPD(NCOUP)=ZRAD
          MRSURF=0
          MTSURF=0
          MASURF=0
          NINCX=0
          NINCZ=0
          IF (NLTRC) THEN
            WRITE (6,*) 'REENTRY FOUND, MPSURF,ZRAD = ',MPSURF,ZRAD
            WRITE (6,*) 'IRCELL ',IRCELL
          ENDIF
          Y00=Y00+VELY*ZRD
          NPCELL=JUPC(NCOUP)
          IPCELL=NPCELL
          GOTO 5000
        ELSE
C  NO REENTRY FOUND
          NCOUP=1
          JUPC(1)=1
          ALPD(1)=ZRAD
          IF (MRSURF.GT.0) THEN
C  CHECK VALID RANGE ON MRSURF
            Y0TEST=Y00+VELY*ZRAD
            IF (NLTRC) WRITE (6,*) 'CHECK VALID RANGE: Y0TEST ',Y0TEST
            IF (Y0TEST.GE.PSURF(1).AND.Y0TEST.LE.PSURF(NP2ND)) THEN
              IPCELL=LEARCA(Y0TEST,PSURF,1,NP2ND,1,'TIMEP 2     ')
            ELSE
              MRSURF=0
              MTSURF=0
              NINCX=0
              NINCZ=0
            ENDIF
          ENDIF
          MPSURF=0
          NINCY=0
          IF (NLTRC) THEN
            WRITE (6,*) 'NO REENTRY FOUND '
          ENDIF
          Y00=Y00+ZRD*VELY
          GOTO 5000
        ENDIF
C
C  PARTICLE IN STANDARD MESH, RADIAL CELL NRCELL
C
2900    CONTINUE
        Y000=Y00+VELY*ZRAD
        IF (NLTRC) WRITE (6,*) 'Y000,Y00,ZRAD ',Y000,Y00,ZRAD
C
        DUM=0.
        NCOUP=1
C
C  J1: CELL INDEX
C  J2: SURFACE INDEX
        J1=NPCELL
        IF (VELY.LT.0.) THEN
          INCY=0
          NINCY=-1
          IF (NLSRFY) J1=MPSURF-1
        ELSE
          INCY=1
          NINCY=1
          IF (NLSRFY) J1=MPSURF
        ENDIF
        J2=J1+INCY
C
        NLSRFY=.FALSE.
C
3000    CONTINUE
        IF (J2.LE.0.OR.J2.GT.NP2ND) THEN
          WRITE (6,*) 'ERROR IN TIMEP ',J2,J1,VELY
c slmod begin (sl)
c...      ADD A QUALIFIER CHECK HERE (GRIDOPT OR SOMETHING):
          WRITE (0,*) 'ERROR 1 IN TIMEP (NOT STOPPING CODE)'
          ZRAD=1.0E+30
          RETURN
c
c          CALL EXIT
c slmod end
        ENDIF
C  TIME FROM Y00 TO PSURF
        IF (MPSURF.EQ.J2) THEN
          J1=J1+NINCY
          J2=J1+INCY
          GOTO 3000
        ENDIF
        DY=(PSURF(J2)-Y00)
        BB=VELY+EPS60
        F=DY/BB
        IF (NLTRC) WRITE (6,*) 'J2,F,DY ',J2,F,DY
        IF (F.LE.ZRAD) THEN
          JUPC(NCOUP)=J1
          LUPC(NCOUP)=J2
          ALPD(NCOUP)=F-DUM
          DUM=F
C  STOP HISTORY AT NON DEFAULT STANDARD SURFACE J2
          ITEST=INMP2I(0,J2,0)
          IN=ITEST+NLIM
!pb          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
          IF ((ityp==3).or.(ITEST.NE.0.AND.ILIIN(IN).NE.0)) THEN
            ZRAD=F
            ISRFCL=0
            NPCELL=J1
            MPSURF=J2
            MRSURF=0
            NINCX=0
            NINCZ=0
            MTSURF=0
            MASURF=0
            IPCELL=NPCELL
            Y00=PSURF(J2)
            GOTO 5000
          ENDIF
C  NEXT CELL
          J1=J1+NINCY
          J2=J1+INCY
          NCOUP=NCOUP+1
          GOTO 3000
        ENDIF
C
C  LAST CELL
C
3100    CONTINUE
        IF (NLTRC) WRITE (6,*) 'LAST CELL ',ZRAD,DUM
        IF (MPSURF.EQ.0) THEN
          NPCELL=J1
        ELSE
          Y0T=Y00+VELY*ZRAD
          NPCELL=LEARCA(Y0T,PSURF,1,NP2ND,1,'TIMEP 3     ')
        ENDIF
        MPSURF=0
        JUPC(NCOUP)=J1
        ALPD(NCOUP)=ZRAD-DUM
        IPCELL=NPCELL
        Y00=Y000
        GOTO 5000
C
      ENDIF
C
5000  CONTINUE
      DO 5100 J=1,NCOUP
        CLPD(IADD+J)=ALPD(J)
        NUPC(IADD+J)=(JUPC(J)-1)+(NTCELL-1)*NP2T3
        NCOUNP(IADD+J)=JUPC(J)
        NCOUNS(IADD+J)=NTCELL
        IF (ALPD(J).LE.0..OR.JUPC(J).LE.0.OR.JUPC(J).GE.NP2ND) THEN
          WRITE (6,*) 'ERROR DETECTED IN TIMEP '
          WRITE (6,*) 'J,IADD+J,ALPD,JUPC ',J,IADD+J,ALPD(J),JUPC(J)
c slmod begin (debug)
          WRITE (6,*) 'ERROR DETECTED IN TIMEP (NP2ND) ',NP2ND
          DO I = 1,NCOUP
            WRITE(6,*) J,NCOUP,JUPC(J),ALPD(J),NP2ND
          ENDDO
c slmod end
        ENDIF
        IF (NLTRC) WRITE (6,*) 'TIMEP ',
     .      J+IADD,CLPD(J+IADD),NUPC(J+IADD),NCOUNP(J+IADD)
5100  CONTINUE
C
      IF (ICOUT.LT.NCOUT.AND.MPSURF.EQ.0) THEN
        ICOUT=ICOUT+1
        IPOLGN=0
        IF (ICOUT.EQ.NCOUT) IPOLGN=IPOLGS
        ZRAD=BLPD(ICOUT)
        if (ncpan.gt.0) then
          jn=0
          do j=ncpan,ncpen
            jn=jn+1
            alpd(jn)=alpd(j)
            jupc(jn)=jupc(j)
            lupc(jn)=lupc(j)
            mupc(jn)=mupc(j)
          end do
        endif
        NTCELL=NCOUNT(ICOUT)
        IADD=IADD+NCOUP
        ncoup=jn
        IF (NLTRC)
     .  WRITE (6,*) 'NEXT TOR. CELL: ZRAD,NTCELL,IADD ',ZRAD,NTCELL,IADD
        GOTO 10000
      ENDIF
C
      NCOU=IADD+NCOUP
C
      SUM=0.
      DO 5110 ICOU=1,NCOU
        SUM=SUM+CLPD(ICOU)
        NCOUNT(ICOU)=NCOUNS(ICOU)
C       WRITE (6,*) 'ICOU,NCOUNT,CLPD ',ICOU,NCOUNT(ICOU),CLPD(ICOU)
5110  CONTINUE
      IF (MPSURF.EQ.0.AND.ABS(SUM-ZRADS).GT.EPS10) THEN
        WRITE (6,*) 'ERROR IN TIMEP: NPANU,SUM,ZRADS ',
     .                               NPANU,SUM,ZRADS
        WRITE (6,*) 'TRY TO KILL PARTICLE ASAP '
        SUM=-1.0D0
      ENDIF
C
      ZRAD=SUM
c slmod begin (debug)
      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(5X,A,5I4,3X,1P,E12.5)')
     .    'TIMEP: Output (MPSURF NP,IP,IR,IPOLGN ZRAD) ',
     .    MPSURF,NPCELL,IPCELL,IRCELL,IPOLGN,ZRAD
        WRITE(6,'(5X,A)') 'TIMEP: '
      ENDIF
c slmod end
      RETURN
C
991   CONTINUE
      WRITE (6,*) 'REENTRANCE FROM VACUUM REGION AND LEVGEO=1 IN TIMEP'
      CALL EXIT
      END
c === ROUTINE: timer
C
C  FULL EIRENE GEOMETRY BLOCK  (GEO3D)
C
C
      SUBROUTINE TIMER (PT)
C
C  THIS SUBROUTINE CALCULATES INTERSECTION TIMES IN THE STANDARD
C  MESH "RSURF" (X- OR RADIAL DIRECTION)
C
C  INPUT:
C       NRCELL = CELL NUMBER FOR WHICH NEXT INTERSECTION
C                TIME IS TO BE CALCULATED (I.E. NOT NECESSARLY THE
C                CELL WHICH CONTAINS THE STARTING POINT X0,Y0,Z0)
C                NRCELL=0 IF PARTICLE OUTSIDE STANDARD MESH
C       NJUMP = 0 MEANS: THIS IS THE FIRST CALL OF TIMER FOR THIS TRACK
C                        IN THIS CASE , NLSRFX MUST BE KNOWN
C           NLSRFX = .TRUE. :PARTICLE ON A SURFACE, IN THIS CASE
C                            THE SUBROUTINE NEEDS 'MRSURF'
C                            MRSURF = NUMBER OF THIS SURFACE
C
C           NLSRFX = .FALSE.:PARTICLE NOT ON A SURFACE
C
C           X0,Y0,Z0 = STARTING POINT OF THIS TRACK
C           VELX,VELY,VELZ = VELOCITY OF PARTICLE
C       NJUMP = 1  X0,Y0,Z0,VELX,VELY,VELZ ARE THE SAME
C                  AS IN THE PREVIOUS CALL
C       NJUMP = 2  ONLY VELX,VELY,VELZ ARE THE SAME
C                  AS IN THE PREVIOUS CALL, I.E. PARTICLE HAS
C                  BEEN MOVED BUT VELOCITY HAS NOT BEEN CHANGED
C                  (TO BE WRITTEN)
C  IF (LEVGEO=1 OR LEVGEO=2):
C       TIMINT NE 0. MEANS: INTERSECTION TIME IS KNOWN FROM AN EARLIER
C                CALL AND ITS VALUE IS TIMINT.
C                OTHERWISE (TIMINT=0) IT HAS TO BE CALCULATED IN THIS
C                CALL
C  IF (LEVGEO=3):
C       TIMINT NE 0. MEANS: INTERSECTION TIME IS KNOWN FROM AN EARLIER
C                CALL AND IS TO BE FOUND IN THE ARRAYS TIMPOL,IIMPOL.
C                OTHERWISE (TIMINT=0) IT HAS TO BE CALCULATED IN THIS
C                CALL
C  IF (LEVGEO=4):
C       TIMINT NE 0. MEANS: NOTHING
C  IF (LEVGEO=5):
C       TIMINT NE 0. MEANS: NOTHING
C  IF (LEVGEO=6):
C       TIMINT NE 0. MEANS: NOTHING
C
C  OUTPUT :
C       NJUMP = 1
C       MRSURF = INDEX OF NEXT SURFACE ALONG TRACK
C              = 0 IF NO NEXT SURFACE IS FOUND
C       PT = TIME TO REACH THIS SURFACE
C          = 1.D30 IF NO NEXT SURFACE IS FOUND
C       TIMINT(MRSURF) NE 0 INDICATES FURTHER INTERSECTION TIMES FOUND
C                      IN THIS CALL, WHICH MAY BE USED IN A LATER CALL
C       NINCX = INDEX FOR DIRECTION IN GRID "RSURF": +1 OR -1
C             = 0 IF NO NEXT SURFACE FOUND
C       NRCELL:  NUMBER OF FINAL RADIAL CELL (IE. NOT MODIFIED)
C       IRCELL: NOT NEEDED ON RADIAL SURFACE. SET IRCELL=NRCELL
C
C  ADDITIONALLY IF NLPLG:
C  INPUT  :
C       IPOLG  = INDEX ON POLYGON OF INITIAL POINT X=(X0,Y0,Z0)
C  OUTPUT :
C       IPOLGN = INDEX ON POLYGON OF THE POINT X+PT*VEL
C         ( IF NO VALID INTERSECTION FOUND IN THIS CALL, THEN
C           IPOLGN=IPOLGO IS RETURNED, THE INDEX OF THE LAST
C           VALID POINT OF INTERSECTION FOUND IN EARLIER CALLS
C           (WHICH MAY BE THE INPUT VALUE IPOLG ITSELF) )
C
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CTSURF
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTETRA
      USE COMPRT
      USE CLGIN
      USE CTRIG

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: PT
      REAL(DP) :: A(3,3), B(3), AB(3,3)
      REAL(DP) :: SIG1, SIG2, DET, XMU, XETA, SARRUS, AX, AY, V, TIMTET,
     .          ZZ, RICHTY, HELP, TM, RICHTX, YY, XX, S3, TIMT, TIM,
     .          S2, A1, A2, A3, B1, B2, B3, S1, V2, V3, ZR, ZEP1, XETAT,
     .          ZSQRT, PS, VELYQ, XX0, VVELX, VELXQ, ZT1, ZT2, YVY,
     .          ZC1, DXA, TST, XA, ZB2, ZAB, ZAB2, ZB, Z0TEST, XMUT,
     .          Y0TEST, ZA, T, PT1, PT2, PT3, PT4, V1, ESURF, DETT,
     .          XTEST, YTEST, X0SURF, DSRF, Y0Q, T1, T2, T3, T4, TPB
      INTEGER :: IRICH(2,4), ITSIDE(3,4)
      INTEGER :: IZELLO, NTIMT, IPOLGOO, IOB, I1, I2, ISW, KAN, KEN,
     .           ILLZ, IHELP, IZELL, NTMS, MXSF, NTMZ, NRMSRF, ICOS,
     .           IERR, IRS, NEWCEL, ITET, IT, IL, IS, NRI, MS, IR,
     .           ICALL, ITFRST, ISTS, MMSURF, ICOUP, J, K, I, JPOL,
     .           MPOL, IPOLGO, LEARC2, IP, ICELLR, MSAVE
      LOGICAL :: LCUT(N2NDPLG)
      LOGICAL :: LNGB1, LNGB2, LNGB3, LNGB4,
     .           LCT1, LCT2, LCT3, LCT4, BITGET
c slmod begin
      INTEGER :: IR1, VP1, ILOOP, AP
      REAL(DP) :: VPLX1,VPLY1,VPLX3,VPLY3,VVTX2,VVTY2,VVTX4,VVTY4
c...  Temp:
      LOGICAL REMARK
      DATA REMARK /.TRUE./
c slmod end
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/
      DATA IRICH / 1, -3,
     .             4,  1,
     .             5,  2,
     .             6,  3 /
      DATA ICALL /0/
      DATA ITFRST /0/
      SAVE
C
      IF (NLTRC) THEN
        CALL LEER(1)
        WRITE (6,*) 'TIMER, INIT.: NRCELL,NLSRFX,MRSURF,TT,TL'
        WRITE (6,*)                NRCELL,NLSRFX,MRSURF,TT,TL
        WRITE (6,*) '              NPCELL,NLSRFY,MPSURF'
        WRITE (6,*)                NPCELL,NLSRFY,MPSURF
      ENDIF
C
c slmod begin (sl)
c...  Count the number of search iterations for radial surface collisions
c     so that infinite loops can be detected:
      ILOOP=0
c slmod end
      IRCELL=NRCELL
      PT=1.D30
      IF (ABS(VELZ).EQ.1.D0) THEN
        IPOLGN=IPOLG
        NCOUP=1
        ALPD(NCOUP)=PT
        LUPC(NCOUP)=0
        MUPC(NCOUP)=0
        JUPC(NCOUP)=IPOLGN
        RETURN
      ENDIF
C
C-------------------------------------------------------------------
      IF (LEVGEO.GT.1) GOTO 100
C
C****SLAB-MODEL IN X DIRECTION
C
      IF (ABS(VELY).EQ.1.D0) THEN
        IPOLGN=IPOLG
        RETURN
      ENDIF
C
      IF (NJUMP.EQ.0) THEN
        XA=X0
        IF (NLSRFX) XA=RSURF(MRSURF)
        NINCX=1
        IF (VELX.LT.0.) NINCX=-1
C
C   RESET BRANCH MARKER
        NJUMP=1
        NLSRFX=.FALSE.
      ENDIF
C
      MRSURF=0
      IF (NRCELL.GT.0) THEN
        IR=NRCELL
        IF (NINCX.EQ.1) IR=IR+1
        DXA=RSURF(IR)-XA
        TST=DXA/(VELX+EPS60)
        IF (NLTOR.AND.NLTRZ) THEN
          Z0TEST=Z0+TST*VELZ
          IF (ZSURF(1).GT.Z0TEST.OR.ZSURF(NT3RD).LT.Z0TEST) GOTO 5
        ENDIF
        IF (NLPOL) THEN
          Y0TEST=Y0+TST*VELY
          IF (PSURF(1).GT.Y0TEST.OR.PSURF(NP2ND).LT.Y0TEST) GOTO 5
        ENDIF
        PT=TST
        MRSURF=IR
5       CONTINUE
      ELSE
C  TRY TO FIND REENTRY SURFACE. CHECK ONLY NON DEFAULT RADIAL SURFACES
        DO 10 IR=1,NR1ST
          IF (INMP1I(IR,0,0).NE.0) THEN
            DXA=RSURF(IR)-XA
            TST=DXA/(VELX+EPS60)
            IF (TST.LE.0..OR.TST.GT.PT) GOTO 10
            IF (NLTOR.AND.NLTRZ) THEN
              Z0TEST=Z0+TST*VELZ
              IF (ZSURF(1).GT.Z0TEST.OR.ZSURF(NT3RD).LT.Z0TEST) GOTO 10
            ENDIF
            IF (NLPOL) THEN
              Y0TEST=Y0+TST*VELY
              IF (PSURF(1).GT.Y0TEST.OR.PSURF(NP2ND).LT.Y0TEST) GOTO 10
            ENDIF
            PT=TST
            MRSURF=IR
          ENDIF
10      CONTINUE
        IF (MRSURF.EQ.0) NINCX=0
      ENDIF
      IF (NLTRC) THEN
        WRITE (6,*) 'TIMER, OUT: PT,MRSURF,NINCX'
        WRITE (6,*)              PT,MRSURF,NINCX
      ENDIF
      RETURN
C
C
100   CONTINUE
C
C---------------------------------------------------------------------
      IF (LEVGEO.GT.2) GOTO 6000
C
      IF (NLELL) GOTO 1000
C
C****CIRCULAR MESH   (CONCENTRIC)
C
      IF (NLTRC) THEN
        WRITE (6,*) 'NJUMP,NINCX,NRCELL '
        WRITE (6,*)  NJUMP,NINCX,NRCELL
      ENDIF
      IF (NJUMP.EQ.0) THEN
        ZA=VELX*VELX+VELY*VELY
        ZB=X0*VELX+Y0*VELY
        ZB2=ZB*ZB
        ZAB=-ZB/(ZA+EPS60)
        ZAB2=ZA/(ZB2+EPS60)
C
C  TEST FOR DIRECTION
        NINCX=1
        IF (ZB.LT.0) NINCX=-NINCX
C
        IF (NLSRFX) THEN
          ZC1=RQ(MRSURF)
C  IN THIS CASE: DON'T TRUST THE VALUE OF NRCELL. RECOMPUTE FROM MRSURF
          IF (NINCX.EQ.1) NRCELL=MRSURF
          IF (NINCX.EQ.-1) NRCELL=MRSURF-1
        ELSE
          ZC1=X0*X0+Y0*Y0
        ENDIF
      ENDIF
C
      IF (NRCELL.GT.0) THEN
C  PARTICLE INSIDE STANDARD MESH
C  FIND NEXT SURFACE MRSURF
        NRI=0
        MRSURF=NRCELL
        IF (NINCX.EQ.1) MRSURF=MRSURF+1
        GOTO 204
      ENDIF
C
C  PARTICLE OUTSIDE STANDARD MESH
C  FIND NEXT SURFACE MRSURF
      NRI=1
      PS=PT
      MS=0
      IS=0
200   DO 205 IR=NRI,NR1ST
        IF (INMP1I(IR,0,0).NE.0) THEN
          MRSURF=IR
          GOTO 204
        ENDIF
205   CONTINUE
      PT=PS
      MRSURF=MS
      RETURN
C
C  CHECK SURFACE: MRSURF
C
204   IF (MRSURF.EQ.1) GO TO 201
      IF (TIMINT(MRSURF).GT.0.0) GO TO 303
203   ZR=RQ(MRSURF)
C  CHECK FOR ROOT
      ZEP1=ZAB2*(ZC1-ZR)
      IF (ZEP1.LT.1.0) GO TO 202
C  NO ROOT - PATH IN OTHER DIRECTION. THIS MUST BE THE
C  OUTWARD DIRECTION NOW, DUE TO CONVEXITY OF THE MESH
      IF (NRI.GT.0) GOTO 310
201   CONTINUE
      NINCX=1
      MRSURF=MRSURF+1
      IF (NJUMP.EQ.0) GO TO 203
      GO TO 303
202   CONTINUE
C
C  RESET BRANCH MARKER
      NJUMP=1
C
C  INTERSECTION TIMES
      ZSQRT=1.0+SQRT(1.0-ZEP1)
      ZT1=ZAB*ZEP1/ZSQRT
      ZT2=ZAB*ZSQRT
      NLSRFX=.FALSE.
C
CL              3.         RETURN RESULT
C
300   CONTINUE
C
C  CHECK SIGN OF RESULT
      IF(ZEP1.GT.0.0) GO TO 301
C  ONE ROOT NEGATIVE OR ZERO
      PT=MAX(ZT1,ZT2)
      GOTO 310
C
301   CONTINUE
      PT=ZT1
      TIMINT(MRSURF)=ZT2
      GOTO 310
C
C  ROOT ALREADY KNOWN
303   PT=TIMINT(MRSURF)
      GOTO 310
C
310   CONTINUE
      IF (NRI.EQ.0) THEN
        RETURN
      ELSE
        IF (PT.GT.0.AND.PT.LT.PS) THEN
          IS=MRSURF-MAX(0,NINCX)
          MS=MRSURF
          PS=PT
        ENDIF
        NRI=IR+1
        GOTO 200
      ENDIF
C
1000  CONTINUE
C
C****ELLIPTICAL MESH  (NOT NECESSARILY CONCENTRIC OR CONFOCAL)
C
      IF (NRCELL.LE.0) THEN
        WRITE (6,*) 'NRCELL.LE.0 IN TIMER: NOT READY '
        CALL EXIT
      ENDIF
C  COEFFICIENTS OF QUADRATIC
      IF (NJUMP.EQ.0) THEN
        YVY=Y0*VELY
        VELXQ=VELX*VELX
        VELYQ=VELY*VELY
        XX0=X0
        VVELX=VELX
        IF (NLSRFX) THEN
          X0SURF=X0-EP1(MRSURF)
          DSRF=ELLQ(MRSURF)
          ZR=RQ(MRSURF)
          Y0Q=(ZR-X0SURF*X0SURF)*DSRF
          MSAVE=MRSURF
C  IN THIS CASE: DON'T TRUST THE VALUE OF NRCELL. RECOMPUTE FROM MRSURF
          ZB=X0SURF*VELX+Y0*VELY/ELLQ(MRSURF)
C  TEST FOR DIRECTION. CASE: NLSRFX, MRSURF KNOWN
          NINCX=1
          IF (ZB.LT.0)     NINCX =-NINCX
          IF (NINCX.EQ.1)  NRCELL=MRSURF
          IF (NINCX.EQ.-1) NRCELL=MRSURF-1
C  NEXT SURFACE
          MRSURF=NRCELL
          IF (NINCX.EQ.1) MRSURF=MRSURF+1
          GOTO 1100
        ELSE
          MSAVE=0
          Y0Q=Y0*Y0
C   TEST FOR DIRECTION. CASE:.NOT.NLSRFX, NRCELL KNOWN
          NINCX=-1
          MRSURF=NRCELL
          ZB=(XX0-EP1(MRSURF))*VVELX+YVY/ELLQ(MRSURF)
          IF(ZB.LE.0.) GOTO 1200
          NINCX=1
          MRSURF=MRSURF+1
          GOTO 1100
        ENDIF
      ENDIF
C
2000  CONTINUE
      MRSURF=NRCELL
      IF (NINCX.EQ.1) MRSURF=MRSURF+1
1100  ZB=(XX0-EP1(MRSURF))*VVELX+YVY/ELLQ(MRSURF)
C
1200  IF(MRSURF.EQ.1) GO TO 2010
      IF(TIMINT(MRSURF).GT.0.0) GO TO 3030
      ESURF=XX0-EP1(MRSURF)
      DSRF=ELLQ(MRSURF)
      ZA=VELXQ+VELYQ/DSRF
      ZB2=ZB*ZB
      ZAB=-ZB/(ZA+EPS60)
      ZAB2=ZA/(ZB2+EPS60)
      ZC1=ESURF*ESURF+Y0Q/DSRF
      ZR=RQ(MRSURF)
      ZEP1=ZAB2*(ZC1-ZR)
C  CHECK FOR ROOT
      IF(ZEP1.LT.1.0) GO TO 2020
C  NO ROOT - PATH IN OTHER DIRECTION
C  MUST BE OUTWARD, BECAUSE OF CONVEXITY OF MESH
2010  CONTINUE
      NINCX=1
      MRSURF=MRSURF+1
      IF (NJUMP.EQ.0.) GO TO 1100
      GO TO 3030
2020  CONTINUE
C
C  RESET BRANCH MARKER
      NJUMP=1
C
C  INTERSECTION TIMES
      ZSQRT=1.0+SQRT(1.0-ZEP1)
      IF (MRSURF.EQ.MSAVE) THEN
        ZT1=0.
        ZEP1=0.
      ELSE
        ZT1=ZAB*ZEP1/ZSQRT
      ENDIF
      ZT2=ZAB*ZSQRT
      NLSRFX=.FALSE.
C
C---------------------------------------------------------------------
CL              3.         RETURN RESULT
C
3000  CONTINUE
C
C  CHECK SIGN OF RESULT
      IF(ZEP1.GT.0.0) GO TO 3010
C  ONE ROOT NEGATIVE OR ZERO
      PT=MAX(ZT1,ZT2)
      GOTO 3100
C
3010  CONTINUE
C  BOTH ROOTS POSITIVE - RETURN ONE AND SAVE OTHER
      PT=ZT1
      TIMINT(MRSURF)=ZT2
      GOTO 3100
C
C  ROOT ALREADY KNOWN
3030  PT=TIMINT(MRSURF)
      GOTO 3100

C
3100  CONTINUE
      XTEST=X0+PT*VELX
      YTEST=Y0+PT*VELY
      IF (NLPOL) IPOLGN=LEARC2(XTEST,YTEST,NRCELL,NPANU,'TIMER')
      IF (NLTRC) WRITE (6,*) 'IPOLGN,NRCELL FROM TIMER ',IPOLGN,NRCELL
      IF (NRI.EQ.0) THEN
        RETURN
      ELSE
        IF (PT.GT.0.AND.PT.LT.PS) THEN
          IS=MRSURF-MAX(0,NINCX)
          MS=MRSURF
          PS=PT
        ENDIF
        NRI=IR+1
        GOTO 2000
      ENDIF
C
C     POLYGONS
C
6000  CONTINUE
C
      IF (LEVGEO.GT.3) GOTO 8000
C
C  POLYGON MESH
C
c slmod begin (debug)
      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(5X,A)')
     .    'TIMER: Input  (NRCELL,NPCELL,NJUMP,NLSRFX,NLSRFY,'//
     .    'NLSRFT MRSURF,MPSURF,IPOLG) '
        WRITE(6,'(5X,3I4,2X,3L2,3I5)')
     .    NRCELL,NPCELL,NJUMP,NLSRFX,NLSRFY,NLSRFT,MRSURF,MPSURF,IPOLG
      ENDIF
c slmod end
      IF (NRCELL.NE.0) THEN
c slmod begin (sl)
c...
        AP=0
c slmod end
        IF (NLTRC) WRITE (6,*) ' TIMER: IN NEIGHBOR PART'
        LNGB1=.TRUE.
        LNGB2=.TRUE.
        LNGB3=.TRUE.
        LNGB4=.TRUE.
C  SET INDICES OF STARTING CELL
        IR=NRCELL
        IP=NPCELL
        IF (NJUMP.EQ.1) THEN
          IR=ICELLR
          IP=IPOLGN
          IF (NINCX.EQ.1) LNGB1=.FALSE.
          IF (NINCX.EQ.-1) LNGB3=.FALSE.
        ENDIF
        IF (NLSRFX) THEN
          NLSRFX=.FALSE.
          IF (NRCELL.EQ.MRSURF) THEN
            LNGB1=.FALSE.
          ELSE
            LNGB3=.FALSE.
          ENDIF
        ENDIF
        IF (NLSRFY) THEN
          IF (IPOLG.EQ.MPSURF) THEN
            IF (NPCELL.EQ.MPSURF) THEN
              LNGB4=.FALSE.
              IP=NGHPLS(2,IR,MPSURF)
            ELSE
              LNGB2=.FALSE.
              IP=NGHPLS(4,IR,MPSURF)
            ENDIF
          ELSE
            LNGB2=.FALSE.
          ENDIF
c slmod begin (sl)
c...BUG? FIX: Not sure what is going on here, but it seems as if
c        the particle intersects the target surface after
c        being launched.  The really strange thing is that 
c        it is still closer to the next poloidal surface than the 
c        target surface, which doesn't make much sense.  However,
c        this only affects a few particles, so I am not
c        spending much time on it here (although, I think something
c        similar may have been happening to an ASDEX grid that Karl
c        was using when David was there, where EIRENE seemed
c        to be launching a particle from behind the target).  Anyway,
c        I am not sure where the problem is, in DIVIMP or EIRENE, but
c        I am putting in this rather specific fix for now:
        ELSEIF (NLSRFT.AND.NPCELL.EQ.1.AND.NRCELL.EQ.11) THEN
          LNGB4=.FALSE.
          WRITE(0,*) 'WARNING: CRUDE FIX TO TARGET CELL PROBLEM',NPANU
c slmod end
        ENDIF
C
        NCOUP=0
C  CALCULATE INTERSECTIONS OF FLIGHT WITH CELL BOUNDARIES
6001    CONTINUE
c slmod begin (sl)
        ILOOP=ILOOP+1
        IF (IP.EQ.0) THEN
          WRITE(0,*) 'ERROR: IP = 0'
          WRITE(6,*) 'ERROR: IP = 0'
          WRITE(6,*) '       NPANU = ',npanu
          PT = -1.0D0
c          PT = 1.D30
          RETURN
        ENDIF
        IF (ILOOP.GT.N2ND) THEN
          WRITE(0,*) 'ERROR: LIKELY AN INFINITE LOOP IN TIMER'
          WRITE(6,*) 'ERROR: LIKELY AN INFINITE LOOP IN TIMER'
          WRITE(6,*) '       NPANU = ',npanu
          PT=-1.0D0
c          PT=1.D30
          RETURN
        ENDIF
c slmod end
c slmod begin
        IF (GRIDOPT.EQ.1) THEN
c  TURN OFF PROCESSING OF VIRTUAL (ZERO VOLUME) CELLS
          IF (RVRTAG(IR,IP).EQ.1) THEN
            LNGB1=.FALSE.
            LNGB3=.FALSE.
          ENDIF
          IF (PVRTAG(IR,IP+1).EQ.1) LNGB2=.FALSE.
          IF (PVRTAG(IR,IP  ).EQ.1) LNGB4=.FALSE.
        ENDIF
c slmod end
        IF (NLTRC) WRITE (6,*) ' IR,IP,LNGB1,LNGB2,LNGB3,LNGB4',
     .                           IR,IP,LNGB1,LNGB2,LNGB3,LNGB4
        T1=-1.D30
        T2=-1.D30
        T3=-1.D30
        T4=-1.D30
        PT1=-1.D30
        PT2=-1.D30
        PT3=-1.D30
        PT4=-1.D30
c slmod begin
        IF (GRIDOPT.EQ.1) THEN
          IF (LNGB1) THEN
c...        VPLX and VPLY are re-computed on the fly here, and elsewhere in the
c           geometry routines.  They should be replaced with arrays
c           that are calculated in the startup routines:
            VPLX1=XVERT(IR,IP+1,1)-XVERT(IR,IP,1)
            VPLY1=YVERT(IR,IP+1,1)-YVERT(IR,IP,1)
            T1=((XVERT(IR,IP,1)-X0)*VELY-
     .          (YVERT(IR,IP,1)-Y0)*VELX)/
     .         (VELX*VPLY1-VELY*VPLX1+EPS60)
          ENDIF
          IF (LNGB3) THEN
            VPLX3=XVERT(IR,IP+1,2)-XVERT(IR,IP,2)
            VPLY3=YVERT(IR,IP+1,2)-YVERT(IR,IP,2)
            T3=((XVERT(IR,IP,2)-X0)*VELY-
     .          (YVERT(IR,IP,2)-Y0)*VELX)/
     .         (VELX*VPLY3-VELY*VPLX3+EPS60)
          ENDIF
          IF (LNGB2)
     .      T2=((XVERT(IR,IP+1,1)-X0)*VELY-
     .          (YVERT(IR,IP+1,1)-Y0)*VELX)/
     .         (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
          IF (LNGB4)
     .      T4=((XVERT(IR,IP,1)-X0)*VELY-
     .          (YVERT(IR,IP,1)-Y0)*VELX)/
     .         (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)
c...      Temp:
          IF (LNGB1) THEN
            DUM=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .         (VELX*VPLY(IR,IP)-VELY*VPLX(IR,IP)+EPS60)
            CALL CHECKNUM('ZAKA1',IR,IP,T1,DUM)
          ENDIF
          IF (LNGB2) THEN
            DUM=((XPOL(IR,IP+1)-X0)*VELY-(YPOL(IR,IP+1)-Y0)*VELX)/
     .          (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
            CALL CHECKNUM('ZAKA2',IR,IP,T2,DUM)
          ENDIF
          IF (LNGB3) THEN
            DUM=((XPOL(IR+1,IP)-X0)*VELY-(YPOL(IR+1,IP)-Y0)*VELX)/
     .          (VELX*VPLY(IR+1,IP)-VELY*VPLX(IR+1,IP)+EPS60)
            CALL CHECKNUM('ZAKA3',IR,IP,T3,DUM)
          ENDIF
          IF (LNGB4) THEN
            DUM=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .          (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)
            CALL CHECKNUM('ZAKA4',IR,IP,T4,DUM)
          ENDIF

        ELSE
          IF (LNGB1)
     .      T1=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .         (VELX*VPLY(IR,IP)-VELY*VPLX(IR,IP)+EPS60)
          IF (LNGB2)
     .      T2=((XPOL(IR,IP+1)-X0)*VELY-(YPOL(IR,IP+1)-Y0)*VELX)/
     .          (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
          IF (LNGB3)
     .      T3=((XPOL(IR+1,IP)-X0)*VELY-(YPOL(IR+1,IP)-Y0)*VELX)/
     .         (VELX*VPLY(IR+1,IP)-VELY*VPLX(IR+1,IP)+EPS60)
          IF (LNGB4)
     .      T4=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .          (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)
        ENDIF
c
c        IF (LNGB1)
c     .  T1=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
c     .     (VELX*VPLY(IR,IP)-VELY*VPLX(IR,IP)+EPS60)
c        IF (LNGB2)
c     .  T2=((XPOL(IR,IP+1)-X0)*VELY-(YPOL(IR,IP+1)-Y0)*VELX)/
c     .      (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
c        IF (LNGB3)
c     .  T3=((XPOL(IR+1,IP)-X0)*VELY-(YPOL(IR+1,IP)-Y0)*VELX)/
c     .     (VELX*VPLY(IR+1,IP)-VELY*VPLX(IR+1,IP)+EPS60)
c        IF (LNGB4)
c     .  T4=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
c     .      (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)
c slmod end

        IF (LNGB2)

        IF (LNGB3)

        IF (LNGB4)

        IF (NLTRC) WRITE (6,*) ' T1,T2,T3,T4 ',T1,T2,T3,T4
        LCT1 = T1.GE.0.D0 .AND. T1.LE.1.D0
        LCT2 = T2.GE.0.D0 .AND. T2.LE.1.D0
        LCT3 = T3.GE.0.D0 .AND. T3.LE.1.D0
        LCT4 = T4.GE.0.D0 .AND. T4.LE.1.D0
C  CALCULATE TIME OF FLIGHT FROM STARTING POINT X=(X0,Y0,Z0)
C  TO THE BOUNDARY OF THE ACTUELL CELL
c slmod begin
        IF (GRIDOPT.EQ.1) THEN
          VPLX1=XVERT(IR,IP+1,1)-XVERT(IR,IP,1)
          VPLY1=YVERT(IR,IP+1,1)-YVERT(IR,IP,1)
          VPLX3=XVERT(IR,IP+1,2)-XVERT(IR,IP,2)
          VPLY3=YVERT(IR,IP+1,2)-YVERT(IR,IP,2)
          VVTX2=VVTX(IR,IP+1)
          VVTY2=VVTY(IR,IP+1)
          VVTX4=VVTX(IR,IP)
          VVTY4=VVTY(IR,IP)
          IF (ABS(VELX).GT.ABS(VELY)) THEN
            IF (LCT1) PT1=(XVERT(IR,IP  ,1)-X0+VPLX1*T1)/VELX
            IF (LCT2) PT2=(XVERT(IR,IP+1,1)-X0+VVTX2*T2)/VELX
            IF (LCT3) PT3=(XVERT(IR,IP  ,2)-X0+VPLX3*T3)/VELX
            IF (LCT4) PT4=(XVERT(IR,IP  ,1)-X0+VVTX4*T4)/VELX
          ELSE
            IF (LCT1) PT1=(YVERT(IR,IP  ,1)-Y0+VPLY1*T1)/VELY
            IF (LCT2) PT2=(YVERT(IR,IP+1,1)-Y0+VVTY2*T2)/VELY
            IF (LCT3) PT3=(YVERT(IR,IP  ,2)-Y0+VPLY3*T3)/VELY
            IF (LCT4) PT4=(YVERT(IR,IP  ,1)-Y0+VVTY4*T4)/VELY
          ENDIF
c...      Temp:
          IF (ABS(VELX).GT.ABS(VELY)) THEN
            IF (LCT1) THEN
              DUM=(XPOL(IR,IP)-X0+VPLX(IR,IP)*T1)/VELX
              CALL CHECKNUM('ZAQA1',IR,IP,PT1,DUM)
            ENDIF
            IF (LCT2) THEN
              DUM=(XPOL(IR,IP+1)-X0+VVTX(IR,IP+1)*T2)/VELX
              CALL CHECKNUM('ZAQA2',IR,IP,PT2,DUM)
            ENDIF
            IF (LCT3) THEN
              DUM=(XPOL(IR+1,IP)-X0+VPLX(IR+1,IP)*T3)/VELX
              CALL CHECKNUM('ZAQA3',IR,IP,PT3,DUM)
            ENDIF
            IF (LCT4) THEN
              DUM=(XPOL(IR,IP)-X0+VVTX(IR,IP)*T4)/VELX
              CALL CHECKNUM('ZAQA4',IR,IP,PT4,DUM)
            ENDIF
          ELSE
            IF (LCT1) THEN
              DUM=(YPOL(IR,IP)-Y0+VPLY(IR,IP)*T1)/VELY
              CALL CHECKNUM('ZAQA5',IR,IP,PT1,DUM)
            ENDIF
            IF (LCT2) THEN
              DUM=(YPOL(IR,IP+1)-Y0+VVTY(IR,IP+1)*T2)/VELY
              CALL CHECKNUM('ZAQA6',IR,IP,PT2,DUM)
            ENDIF
            IF (LCT3) THEN
              DUM=(YPOL(IR+1,IP)-Y0+VPLY(IR+1,IP)*T3)/VELY
              CALL CHECKNUM('ZAQA7',IR,IP,PT3,DUM)
            ENDIF
            IF (LCT4) THEN
              DUM=(YPOL(IR,IP)-Y0+VVTY(IR,IP)*T4)/VELY
              CALL CHECKNUM('ZAQA8',IR,IP,PT4,DUM)
            ENDIF
          ENDIF

        ELSE
          IF (ABS(VELX).GT.ABS(VELY)) THEN
            IF (LCT1) PT1=(XPOL(IR,IP)-X0+VPLX(IR,IP)*T1)/VELX
            IF (LCT2) PT2=(XPOL(IR,IP+1)-X0+VVTX(IR,IP+1)*T2)/VELX
            IF (LCT3) PT3=(XPOL(IR+1,IP)-X0+VPLX(IR+1,IP)*T3)/VELX
            IF (LCT4) PT4=(XPOL(IR,IP)-X0+VVTX(IR,IP)*T4)/VELX
          ELSE
            IF (LCT1) PT1=(YPOL(IR,IP)-Y0+VPLY(IR,IP)*T1)/VELY
            IF (LCT2) PT2=(YPOL(IR,IP+1)-Y0+VVTY(IR,IP+1)*T2)/VELY
            IF (LCT3) PT3=(YPOL(IR+1,IP)-Y0+VPLY(IR+1,IP)*T3)/VELY
            IF (LCT4) PT4=(YPOL(IR,IP)-Y0+VVTY(IR,IP)*T4)/VELY
          ENDIF
        ENDIF
c
c        IF (ABS(VELX).GT.ABS(VELY)) THEN
c          IF (LCT1) PT1=(XPOL(IR,IP)-X0+VPLX(IR,IP)*T1)/VELX
c          IF (LCT2) PT2=(XPOL(IR,IP+1)-X0+VVTX(IR,IP+1)*T2)/VELX
c          IF (LCT3) PT3=(XPOL(IR+1,IP)-X0+VPLX(IR+1,IP)*T3)/VELX
c          IF (LCT4) PT4=(XPOL(IR,IP)-X0+VVTX(IR,IP)*T4)/VELX
c        ELSE
c          IF (LCT1) PT1=(YPOL(IR,IP)-Y0+VPLY(IR,IP)*T1)/VELY
c          IF (LCT2) PT2=(YPOL(IR,IP+1)-Y0+VVTY(IR,IP+1)*T2)/VELY
c          IF (LCT3) PT3=(YPOL(IR+1,IP)-Y0+VPLY(IR+1,IP)*T3)/VELY
c          IF (LCT4) PT4=(YPOL(IR,IP)-Y0+VVTY(IR,IP)*T4)/VELY
c        ENDIF
c slmod end
        LCT1 = LCT1 .AND. PT1.GE.0.D0
        LCT2 = LCT2 .AND. PT2.GE.0.D0
        LCT3 = LCT3 .AND. PT3.GE.0.D0
        LCT4 = LCT4 .AND. PT4.GE.0.D0
c slmod begin (debug)
        IF (.NOT.(LCT1.OR.LCT2.OR.LCT3.OR.LCT4))THEN
          IF (T1.GE.-1.0D-05.AND.T1.LE.0.0D0) THEN
            WRITE(0,*) 'ALMOST 1: NPANU =',npanu
            WRITE(6,*) 'ALMOST 1: NPANU =',npanu
          ENDIF
          IF (T2.GE.-1.0D-05.AND.T2.LE.0.0D0) THEN
            WRITE(0,*) 'ALMOST 2: NPANU =',npanu
            WRITE(6,*) 'ALMOST 2: NPANU =',npanu
          ENDIF
          IF (T3.GE.-1.0D-05.AND.T3.LE.0.0D0) THEN
            WRITE(0,*) 'ALMOST 3: NPANU =',npanu
            WRITE(6,*) 'ALMOST 3: NPANU =',npanu
          ENDIF
          IF (T4.GE.-1.0D-05.AND.T4.LE.0.0D0) THEN
            WRITE(0,*) 'ALMOST 4: NPANU =',npanu
            WRITE(6,*) 'ALMOST 4: NPANU =',npanu
          ENDIF
        ENDIF

        IF (PRINTOPT.EQ.1) THEN
          WRITE(6,'(5X,A,2I4,I2,1X,4L2,1X,4L2,1X,1P,4E15.7)')
     .      'TIMER: ',
     .      ir,ip,rvrtag(ir,ip),lngb1,lngb2,lngb3,lngb4,
     .      lct1,lct2,lct3,lct4,t1,t2,t3,t4
        ENDIF
c slmod end
        IF (NLTRC) WRITE (6,*) ' LCT1,LCT2,LCT3,LCT4 ',
     .                           LCT1,LCT2,LCT3,LCT4
        IF (NLTRC) WRITE (6,*) ' PT1,PT2,PT3,PT4 ',PT1,PT2,PT3,PT4
        T=MAX(PT1,PT2,PT3,PT4)
        IF (NLTRC) WRITE (6,*) ' T = ',T
C  IF INTERSECTION WITH POLOIDAL BOUNDARY CONTINUE WITH NEIGHBORING CELL
c slmod begin (sl)
c... Not sure what this is for:
          IF (LCT2.AND.LCT4.AND.LCT1) THEN
            WRITE(0,*) 'DISCARDING INTERSECTIONS A'
            WRITE(6,*) 'DISCARDING INTERSECTIONS A'
            LCT2 = .FALSE.
            LCT4 = .FALSE.
            T=MAX(PT1,PT3)
          ENDIF
          IF (LCT2.AND.LCT4.AND.LCT3) THEN
            WRITE(0,*) 'DISCARDING INTERSECTIONS B'
            WRITE(6,*) 'DISCARDING INTERSECTIONS B'
            LCT2 = .FALSE.
            LCT4 = .FALSE.
            T=MAX(PT1,PT3)
          ENDIF
c slmod end
        IF (LCT2.OR.LCT4) THEN
c slmod begin (sl)
          IF (LCT2.AND.LCT4) THEN
c  INTERSECTIONS ARE FOUND FOR BOTH OPPOSING SIDES OF A CELL, KILL
c  THE PARTICLE (THIS CAN HAPPEN WITH DISTORTED CELLS WHERE AN INTERIOR
c  ANGLE IS GREATER THAN 180 DEGREES)
            WRITE(0,*) 'TIMER: 2 INTERSECTIONS, KILLING PARTICLE'
            WRITE(6,*) 'TIMER: 2 INTERSECTIONS, KILLING PARTICLE'
            WRITE(6,*) '       NPANU = ',NPANU
            PT=-1.0D0
c            PT=1.D30
            RETURN
          ENDIF
c slmod end
          LNGB1=.TRUE.
          LNGB2=.TRUE.
          LNGB3=.TRUE.
          LNGB4=.TRUE.
          NCOUP=NCOUP+1
          ALPD(NCOUP)=T
          IF (LCT2) THEN
            LUPC(NCOUP)=IP+1
            MUPC(NCOUP)=1
            JUPC(NCOUP)=IP
            IP=NGHPOL(2,IR,IP)
            LNGB4=.FALSE.
          ENDIF
          IF (LCT4) THEN
            LUPC(NCOUP)=IP
            MUPC(NCOUP)=-1
            JUPC(NCOUP)=IP
            IP=NGHPOL(4,IR,IP)
            LNGB2=.FALSE.
          ENDIF
          ISTS=INMP2I(IR,LUPC(NCOUP),0)
!pb          IF ((.not.NLPOL.or.ISTS.eq.0).and.ip.ne.0) goto 6001
          IF (ityp.ne.3.and.(.not.NLPOL.or.ISTS.eq.0).and.ip.ne.0) 
     .       goto 6001
C  NO NEIGHBORING CELL: PARTICLE HAS HIT A POLOIDAL BOUNDARY OF THE MESH
c slmod begin (debug)
          IF (PRINTOPT.EQ.1)
     .      WRITE(6,'(5X,A)') 'TIMER: PARTICLE HAS HIT A POLOIDAL '//
     .                        'BOUNDARY OF THE MESH'
c slmod end
          MRSURF=0
          PT=1.D30
          NINCX=0
        ELSEIF (LCT1.OR.LCT3) THEN
C  INTERSECTION WITH RADIAL CELL BOUNDARY FOUND
c slmod begin (debug)
          IF (PRINTOPT.EQ.1)
     .      WRITE(6,'(5X,A)') 'TIMER: INTERSECTION WITH RADIAL '//
     .                        'CELL BOUNDARY FOUND'
c slmod end
          PT=T
          IPOLGN=IP
          LNGB1=.TRUE.
          LNGB2=.TRUE.
          LNGB3=.TRUE.
          LNGB4=.TRUE.
          NCOUP=NCOUP+1
          ALPD(NCOUP)=T
          LUPC(NCOUP)=0
          MUPC(NCOUP)=0
          JUPC(NCOUP)=IP
          IF (LCT1) THEN
            MRSURF=IR
            NINCX=-1
            ICELLR=IR-1
          ELSE
            MRSURF=IR+1
            ICELLR=IR+1
            NINCX=1
          ENDIF
        ELSE
C  NO INTERSECTION FOUND
c slmod begin
          IF (GRIDOPT.EQ.1.AND.
     .        PVRTAG(IR,IP).EQ.1.OR.PVRTAG(IR,IP+1).EQ.1) THEN
C  NEED TO ACCOUNT FOR VIRTUAL (ZERO VOLUME) CELLS
            LNGB1=.TRUE.
            LNGB2=.TRUE.
            LNGB3=.TRUE.
            LNGB4=.TRUE.
            IF (NLSRFY.AND.IP.EQ.NPOINT(2,1)-1) THEN
              IF (MPSURF.EQ.NPOINT(2,1)) THEN
                IP=NGHPOL(2,IR,IP)
                LNGB4=.FALSE.
                AP=+1
              ELSE
                IP=NGHPOL(4,IR,IP)
                LNGB2=.FALSE.
                AP=-1
              ENDIF
            ELSEIF (NLSRFY.AND.IP.EQ.NPOINT(1,NPPLG)) THEN
              IF (MPSURF.EQ.NPOINT(1,NPPLG)) THEN
                IP=NGHPOL(4,IR,IP)
                LNGB2=.FALSE.
                AP=-1
              ELSE
                IP=NGHPOL(2,IR,IP)
                LNGB4=.FALSE.
                AP=+1
              ENDIF
            ELSEIF (PVRTAG(IR,IP).EQ.0) THEN
              IP=NGHPOL(2,IR,IP)
              AP=+1
            ELSEIF (PVRTAG(IR,IP+1).EQ.0) THEN
              IP=NGHPOL(4,IR,IP)
              AP=-1
            ELSE
              IP=IP+AP
            ENDIF
            NJUMP=1
            GOTO 6001
          ENDIF
c...      Debug (remove):
          IF (PRINTOPT.EQ.1)
     .      WRITE(6,'(5X,A)') 'TIMER: NO INTERSECTION FOUND '
c slmod end
          IF (NLSRFY.AND.NJUMP.EQ.0) THEN
C  PLAY SAVE: TRY ONCE AGAIN, IF PARTICLE ON POL. SURFACE
            WRITE (6,*) ' NO INTERSECTION IN TIMER. TRY ONCE AGAIN '
            MMSURF=MSURF
            IF (MSURF.GT.NLIM) MMSURF=-MSURF+NLIM
            WRITE (6,*) 'NPANU, MSURF ',NPANU,MMSURF
            IF (.NOT.LNGB4) THEN
              IP=NGHPLS(4,IR,MPSURF)
              LNGB4=.TRUE.
              LNGB2=.FALSE.
            ELSEIF (.NOT.LNGB2) THEN
              IP=NGHPLS(2,IR,MPSURF)
              LNGB2=.TRUE.
              LNGB4=.FALSE.
            ENDIF
            LNGB1=.TRUE.
            LNGB3=.TRUE.
            NJUMP=1
            GOTO 6001
          ENDIF
          WRITE (6,*) ' ERROR: NO INTERSECTION FOUND IN TIMER '
          WRITE (6,*) ' NPANU: ',NPANU
          MRSURF=0
          PT=1.D30
          NINCX=0
          ICELLR=0
          IPOLGN=IP
        ENDIF
        NJUMP=1
        IF (NLTRC) THEN
          WRITE (6,*) ' PT,MRSURF,NINCX,IPOLGN ',
     .                  PT,MRSURF,NINCX,IPOLGN
          WRITE (6,*) 'NCOUP ',NCOUP
          DO ICOUP=1,NCOUP
            WRITE (6,*) 'ICOUP,ALPD(ICOUP) ',ICOUP,ALPD(ICOUP)
          ENDDO
        ENDIF
c slmod begin (debug)
        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(5X,A,2I3,2I3,A,E12.5,2I3)')
     .      'TIMER: Output (NR,NJUMP,1X,MR,IPOLGN,PT,NINCX,IR)',
     .      NRCELL,NJUMP,MRSURF,IPOLGN,'*',PT,NINCX,IRCELL
          WRITE(6,'(5X,A)') 'TIMER: '
        ENDIF
c slmod end
        RETURN
      ENDIF
C
C  PARTICLE OUTSIDE STANDARD MESH, IN ADDITIONAL CELL NACELL
C  NRCELL=0
C
      IF (NLSRFX) THEN
        NLSRFX=.FALSE.
        JPOL=IPOLG
        MPOL=MRSURF
      ELSE
        JPOL=0
        MPOL=0
      ENDIF
C
      IF (NJUMP.EQ.0) IPOLGO=IPOLG
C
      PT=1.D30
      MRSURF=0
      IPOLGN=IPOLGO
C
c slmod begin (debug)
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(5X,A)') 'TIMER: CHECKING FOR RE-ENTRY'
c slmod end
      DO 6100 I=1,NR1ST
        ISTS=INMP1I(I,0,0)
C  SURFACE I IS NOT A NON DEFAULT RADIAL STANDARD SURFACE
        IF (ISTS.EQ.0) GOTO 6100
        IF (TIMINT(I).EQ.0) THEN
          TIMINT(I)=1
          NTIM(I)=0
c slmod begin
          IF (GRIDOPT.EQ.1) THEN
c...        As in LEARC1, the DIVSUR check needs to be replaced with
c           a parameter in the input file that indicates which side of 
c           the ring the surface is on:
            IF (I.GT.DIVSUR) THEN
              IR1=I-1
              VP1=2
            ELSE
              IR1=I
              VP1=1
            ENDIF
          ENDIF
c slmod end
C
C   SEARCH FOR ALL POSSIBLE INTERSECTIONS WITHIN THE CELL
C
          DO 6011 J=1,NRPLG
6011        LCUT(J)=.FALSE.
C
          DO 6012 J=1,NPPLG
            DO 6012 K=NPOINT(1,J),NPOINT(2,J)-1
c slmod begin (sl)
c...note: hope this is okay!
              IF (REMARK.AND.INMP1I(I,K,0).EQ.0) THEN
                 REMARK = .FALSE.
                 WRITE(0,*) '*** QUESTIONABLE USE OF INMP1I ***'
              ENDIF
              IF (GRIDOPT.EQ.0.AND.INMP1I(I,K,0).EQ.0) GOTO 6012
c slmod end
c slmod begin
              IF (GRIDOPT.EQ.1) THEN
                IF (INMP1I(I,K,0).EQ.0) CYCLE
                V1=(YVERT(IR1,K+1,VP1)-Y0)*VELX-
     .             (XVERT(IR1,K+1,VP1)-X0)*VELY
                V2=(YVERT(IR1,K  ,VP1)-Y0)*VELX-
     .             (XVERT(IR1,K  ,VP1)-X0)*VELY
c...            Temp:
                DUM=(YPOL(I,K+1)-Y0)*VELX-(XPOL(I,K+1)-X0)*VELY
                CALL CHECKNUM('ZEWA1',IR,IP,V1,DUM)
                DUM=(YPOL(I,K)-Y0)*VELX-(XPOL(I,K)-X0)*VELY
                CALL CHECKNUM('ZEWA2',IR,IP,V2,DUM)
              ELSE
                V1=(YPOL(I,K+1)-Y0)*VELX-(XPOL(I,K+1)-X0)*VELY
                V2=(YPOL(I,K)-Y0)*VELX-(XPOL(I,K)-X0)*VELY
              ENDIF
c
c              V1=(YPOL(I,K+1)-Y0)*VELX-(XPOL(I,K+1)-X0)*VELY
c              V2=(YPOL(I,K)-Y0)*VELX-(XPOL(I,K)-X0)*VELY
c slmod end
              LCUT(K)=V1*V2.LE.0.
6012      CONTINUE
C
          IF (I.EQ.MPOL) LCUT(JPOL)=.FALSE.
          KAN=ILLZ(NRPLG,LCUT,1)+1
          KEN=NRPLG-ILLZ(NRPLG,LCUT,-1)
C
C   COMPUTE THE FLIGHT TIMES TO THE INTERSECTION POINTS
C
          DO 6013 K=KAN,KEN
            IF (LCUT(K)) THEN
c slmod begin
              IF (GRIDOPT.EQ.1) THEN
                VPLX1=XVERT(IR1,K+1,VP1)-XVERT(IR1,K,VP1)
                VPLY1=YVERT(IR1,K+1,VP1)-YVERT(IR1,K,VP1)
                T1=((XVERT(IR1,K,VP1)-X0)*VPLY1-
     .              (YVERT(IR1,K,VP1)-Y0)*VPLX1)/
     .             (VELX*VPLY1-VELY*VPLX1+EPS60)
c...            Temp:
                DUM=((XPOL(I,K)-X0)*VPLY(I,K)-(YPOL(I,K)-Y0)*VPLX(I,K))
     .              /(VELX*VPLY(I,K)-VELY*VPLX(I,K)+EPS60)
                CALL CHECKNUM('ZEWA4',I,K,T1,DUM)
              ELSE
                T1=((XPOL(I,K)-X0)*VPLY(I,K)-(YPOL(I,K)-Y0)*VPLX(I,K))
     .             /(VELX*VPLY(I,K)-VELY*VPLX(I,K)+EPS60)
              ENDIF
c
c              T1=((XPOL(I,K)-X0)*VPLY(I,K)-(YPOL(I,K)-Y0)*VPLX(I,K))
c     .           /(VELX*VPLY(I,K)-VELY*VPLX(I,K)+EPS60)
c slmod end
              IF (T1.GT.0.) THEN
                NTIM(I)=NTIM(I)+1
                TIMPOL(I,NTIM(I))=T1
                IIMPOL(I,NTIM(I))=K
              ENDIF
            ENDIF
6013      CONTINUE
C
6015      ISW=0
          DO 6020 J=1,NTIM(I)-1
            IF (TIMPOL(I,J).LT.TIMPOL(I,J+1)) THEN
              ISW=ISW+1
              HELP=TIMPOL(I,J)
              TIMPOL(I,J)=TIMPOL(I,J+1)
              TIMPOL(I,J+1)=HELP
              IHELP=IIMPOL(I,J)
              IIMPOL(I,J)=IIMPOL(I,J+1)
              IIMPOL(I,J+1)=IHELP
            ENDIF
6020      CONTINUE
          IF (ISW.GT.0.AND.NTIM(I).GT.2) GOTO 6015
          IF (NLTRC) THEN
            WRITE (6,*) ' SURFACE NO. = ',I
            WRITE (6,*) ' TIMPOL ',(TIMPOL(I,J),J=1,NTIM(I))
            WRITE (6,*) ' IIMPOL ',(IIMPOL(I,J),J=1,NTIM(I))
          ENDIF
C
        ENDIF
C
C  FIND PT, MRSURF, IPOLGN
        IF (NTIM(I).GT.0) THEN
          IF (TIMPOL(I,NTIM(I)).LT.PT) THEN
            MRSURF=I
            PT=TIMPOL(I,NTIM(I))
            IPOLGN=IIMPOL(I,NTIM(I))
          ENDIF
        ENDIF
6100  CONTINUE
C
C  FIND NINCX
      NINCX=0
      IF (MRSURF.NE.0) THEN
        NTIM(MRSURF)=NTIM(MRSURF)-1
        NINCX=SIGN(1._DP,VELX*PLNX(MRSURF,IPOLGN)+
     .        VELY*PLNY(MRSURF,IPOLGN))
      ENDIF
C
      NJUMP=1
      IPOLGO=IPOLGN
C
      IF (NLTRC) WRITE (6,*) 'PT,MRSURF,IPOLGN,NINCX ',
     .                        PT,MRSURF,IPOLGN,NINCX
      RETURN
C
C
8000  CONTINUE
C
      IF (LEVGEO.GT.4) GOTO 10000
C
C   FINITE ELEMENT DISCRETISATION
C
      IF (NLTRC) THEN
        WRITE (6,*) ' TIMER,NJUMP,NRCELL ',NJUMP,NRCELL
        WRITE (6,*) '       NLSRFX,IPOLG ',NLSRFX,IPOLG
      ENDIF
C
      IF (NJUMP.EQ.0) THEN
8001    CONTINUE
        IZELL = NRCELL
        IPOLGO=0
        IF (NLSRFX) THEN
          IPOLGO=IPOLG
        ENDIF
C     ELSEIF (NJUMP.EQ.1) THEN
      ENDIF
        IF (NLTRC) THEN
          WRITE (6,*) ' IZELL,IPOLGO ',IZELL,IPOLGO
          WRITE (6,*) ' X0,Y0 ',X0,Y0
          WRITE (6,*) ' VELX,VELY ',VELX,VELY
          CALL LEER(1)
        ENDIF
        XX = X0
        YY = Y0
        TM=0.
        IF (IZELL.EQ.0) GOTO 9999
C
8020    CONTINUE
        DO 8010,J=1,3
          IF (J.NE.IPOLGO) THEN
            RICHTX = VTRIX(J,IZELL)
            RICHTY = VTRIY(J,IZELL)
            AX = XTRIAN(NECKE(J,IZELL))-XX
            AY = YTRIAN(NECKE(J,IZELL))-YY
            V  = (AX*VELY-AY*VELX)/(RICHTY*VELX-RICHTX*VELY+EPS60)
            IF (V.GE.0..AND.V.LE.1.) THEN
              IF (ABS(VELX).GT.ABS(VELY)) THEN
                T=(AX+V*RICHTX)/VELX
              ELSE
                T=(AY+V*RICHTY)/VELY
              ENDIF
              IF (T .GT. 0.) THEN
                XX = XX + T * VELX
                YY = YY + T * VELY
                TM = TM + T
                TIMINT(IZELL) = TM
C  NUMMER DER GESCHNITTENEN SEITE DES ALTEN DREIECKS
                NTIM(IZELL) = J
                IF (NLTRC) THEN
                  WRITE (6,*) 'IZELL,XX,YY,J ',IZELL,XX,YY,J
                ENDIF
C  ZELLENNUMMER DES NEUEN DREIECKS
                IZELLO=IZELL
                IZELL = NCHBAR(J,IZELLO)
                IF (IZELL .EQ. 0) GOTO 8050
C  SEITENNUMMER DES NEUEN DREIECKS
                IPOLGO = NSEITE(J,IZELLO)
                IF (NLTRC) THEN
                  WRITE(6,*) 'GEHE IN ZELLE ',IZELL,' TM= ',TM
                  WRITE(6,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
                ENDIF
                GOTO 8050
              ENDIF
            ENDIF
          ENDIF
8010    CONTINUE
C
        IF (NLTRC) WRITE (6,*) ' NO INTERSECTION FOUND IN TRIANGLE '
        PT=1.D30
        NINCX=0
C  DURCHSUCHE NACHBARDREIECK
        IF (NLSRFX) THEN
          IZELL=NCHBAR(IPOLG,NRCELL)
          IPOLGO=NSEITE(IPOLG,NRCELL)
          NRCELL=IZELL
          NLSRFX=.FALSE.
          GOTO 8020
        ENDIF
C  KEIN SCHNITTPUNKT GEFUNDEN UND NLSRFX=.FALSE.
        RETURN
C
C
C   THIS IS DONE FOR NJUMP=0 AND NJUMP=1
C
8050  CONTINUE
      NLSRFX=.FALSE.
      NJUMP=1
      MRSURF=NRCELL
      PT=TIMINT(MRSURF)
      TIMINT(MRSURF)=0.
      IPOLGN=NTIM(NRCELL)
      ISTS=ABS(INMTI(IPOLGN,NRCELL))
      IF (ISTS.EQ.0) THEN
        NINCX=NCHBAR(IPOLGN,NRCELL)-NRCELL
      ELSEIF (ISTS.GT.0.AND.
     .        ISTS.LE.NLIM+NSTSI) THEN
C  ON NON DEFAULT SURFACE (ADD. OR STD.) ISTS=INMTI(IPOLGN,NRCELL)
        NINCX=SIGN(1,INMTI(IPOLGN,NRCELL))
      ELSE
        GOTO 9999
      ENDIF
      IF (NLTRC) WRITE (6,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
     .                         NRCELL,MRSURF,PT,NINCX,IPOLGN

      RETURN
C
10000 CONTINUE

      IF (LEVGEO > 5) GOTO 12000
C
C     TETRAHEDRONS
C
      TIMTET = 1.E30
      NTIMT = 0
      IF (NRCELL .NE. 0) THEN
        IF (NJUMP.EQ.0) THEN
          IZELL = NRCELL
          IPOLGO=0
          IF (NLSRFX) THEN
            IPOLGO=IPOLG
          ENDIF
C       ELSEIF (NJUMP.EQ.1) THEN
        ENDIF
        IF (NLTRC) THEN
          WRITE (6,*) ' IZELL,IPOLGO ',IZELL,IPOLGO
          WRITE (6,*) ' X0,Y0,Z0 ',X0,Y0,Z0
          WRITE (6,*) ' VELX,VELY,VELZ ',VELX,VELY,VELZ
        ENDIF

        XX = X0
        YY = Y0
        ZZ = Z0
        TM=0.
        IF (IZELL.EQ.0) THEN
          WRITE (6,*) ' IZELL = 0 IN TIMER '
          CALL EXIT

        END IF
C
C  CHECK TETRAHEDRON IZELL FOR INTERSECTION WITH TRAJECTORY
C
11020   CONTINUE
        DO J=1,4
          IF (J.NE.IPOLGO) THEN
            SIG1=SIGN(1,IRICH(1,J))
            I1=ABS(IRICH(1,J))
            SIG2=SIGN(1,IRICH(2,J))
            I2=ABS(IRICH(2,J))
            A(1:3,1) = (/ VTETX(I1,IZELL), VTETY(I1,IZELL),
     .                    VTETZ(I1,IZELL)/) * SIG1
            A(1:3,2) = (/ VTETX(I2,IZELL), VTETY(I2,IZELL),
     .                    VTETZ(I2,IZELL)/) * SIG2
            A(1:3,3) = (/ -VELX, -VELY, -VELZ /)
            B(1:3) = (/ XX-XTETRA(NTECK(ITSIDE(1,J),IZELL)),
     .                  YY-YTETRA(NTECK(ITSIDE(1,J),IZELL)),
     .                  ZZ-ZTETRA(NTECK(ITSIDE(1,J),IZELL)) /)
cpb            DETT = SARRUS(A)
            DET = A(1,1) * A(2,2) * A(3,3) 
     .          + A(1,2) * A(2,3) * A(3,1)
     .          + A(1,3) * A(2,1) * A(3,2)
     .          - A(3,1) * A(2,2) * A(1,3)
     .          - A(3,2) * A(2,3) * A(1,1)
     .          - A(3,3) * A(2,1) * A(1,2)
            IF (ABS(DET) < 1.D-5) CYCLE
            AB = A
            AB(:,1) = B
cpb            XMUT = SARRUS(AB)/DET
            XMU = ( AB(1,1) * AB(2,2) * AB(3,3) 
     .            + AB(1,2) * AB(2,3) * AB(3,1)
     .            + AB(1,3) * AB(2,1) * AB(3,2)
     .            - AB(3,1) * AB(2,2) * AB(1,3)
     .            - AB(3,2) * AB(2,3) * AB(1,1)
     .            - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
            IF (XMU .GE.0.D0 .AND. XMU .LE.1.D0) THEN
              AB = A
              AB(:,2) = B
cpb              XETAT = SARRUS(AB)/DET
              XETA = ( AB(1,1) * AB(2,2) * AB(3,3) 
     .               + AB(1,2) * AB(2,3) * AB(3,1)
     .               + AB(1,3) * AB(2,1) * AB(3,2)
     .               - AB(3,1) * AB(2,2) * AB(1,3)
     .               - AB(3,2) * AB(2,3) * AB(1,1)
     .               - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
              IF ((XETA.GE.0.D0 .AND. XETA.LE.1.D0) .AND.
     .            (XMU+XETA <= 1.D0)) THEN
                AB = A
                AB(:,3) = B
cpb                TPB = SARRUS(AB)/DET
                T = ( AB(1,1) * AB(2,2) * AB(3,3) 
     .              + AB(1,2) * AB(2,3) * AB(3,1)
     .              + AB(1,3) * AB(2,1) * AB(3,2)
     .              - AB(3,1) * AB(2,2) * AB(1,3)
     .              - AB(3,2) * AB(2,3) * AB(1,1)
     .              - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
                IF (T .GT. 0.) THEN
!            IF ((XMU .GE.0.D0 .AND. XMU .LE.1.D0) .AND.
!     .          (XETA.GE.0.D0 .AND. XETA.LE.1.D0) .AND.
!     .          (XMU+XETA <= 1.D0) .AND. (T .GT. 0.)) THEN
                  XX = XX + T * VELX
                  YY = YY + T * VELY
                  ZZ = ZZ + T * VELZ
                  if (nltrc) then
                     write(6,*) ' i1, i2 ',i1,i2
                     write(6,'(1x,a,3es14.7)') ' velxyz ',velx,vely,velz
                     write(6,'(1x,a,3es14.7)') ' a11, a12, a13 ',a(1,:)
                     write(6,'(1x,a,3es14.7)') ' a21, a22, a23 ',a(2,:)
                     write(6,'(1x,a,3es14.7)') ' a31, a32, a33 ',a(3,:)
                     write(6,'(1x,a,3es14.7)') ' b1, b2, b3 ',b(:)
                     write (6,'(1x,a,3es14.7)') ' det, xmu ',det,xmu
                     write (6,'(1x,a,3es14.7)') ' dett, xmut ',dett,xmut
                     write (6,'(1x,a,3es14.7)') ' xeta, t ',xeta,t
                     write (6,'(1x,a,3es14.7)') ' xetat, tpb ',xetat,tpb
                  end IF
                  TM = TM + T
                  TIMTET = TM
C  NUMMER DER GESCHNITTENEN SEITE DES ALTEN DREIECKS
                  NTIMT = J
                  IF (NLTRC) THEN
                    WRITE (6,*) 'IZELL,XX,YY,ZZ,J ',IZELL,XX,YY,ZZ,J
                  ENDIF
C  ZELLENNUMMER DES NEUEN DREIECKS
                  IZELLO=IZELL
                  IPOLGOO=IPOLGO
                  IZELL = NTBAR(J,IZELLO)
                  PT = TM
                  IF (IZELL .EQ. 0) GOTO 11050
C  SEITENNUMMER DES NEUEN DREIECKS
                  IPOLGO = NTSEITE(J,IZELLO)
                  IF (NLTRC) THEN
                    WRITE(6,*) 'GEHE IN ZELLE ',IZELL,' TM= ',TM
                    WRITE(6,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
                  ENDIF
                  GOTO 11050
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        END DO
C
        IF (NLTRC) WRITE (6,*) ' NO INTERSECTION FOUND IN TETRAHEDRON '
        PT=1.D30
        NINCX=0
C  DURCHSUCHE NACHBARDREIECK
        IF (NLSRFX) THEN
          IZELL=NTBAR(IPOLG,NRCELL)
          IPOLGO=NTSEITE(IPOLG,NRCELL)
          NRCELL=IZELL
          NLSRFX=.FALSE.
          GOTO 11020
        ENDIF
C  KEIN SCHNITTPUNKT GEFUNDEN UND NLSRFX=.FALSE.
        RETURN
C
C
C   THIS IS DONE FOR NJUMP=0 AND NJUMP=1
C
11050   CONTINUE
        NLSRFX=.FALSE.
        NJUMP=1
        MRSURF=NRCELL
        PT=TIMTET
        IPOLGN=NTIMT
        ISTS=ABS(INMTIT(IPOLGN,NRCELL))
        IF (ISTS.EQ.0) THEN
          NINCX=NTBAR(IPOLGN,NRCELL)-NRCELL
        ELSEIF (ISTS.GT.0.AND.
     .          ISTS.LE.NLIM+NSTSI) THEN
C  ON NON DEFAULT SURFACE (ADD. OR STD.) ISTS=INMTI(IPOLGN,NRCELL)
          NINCX=SIGN(1,INMTIT(IPOLGN,NRCELL))
        ELSE
          GOTO 9999
        ENDIF
        IF (NLTRC) WRITE (6,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
     .                           NRCELL,MRSURF,PT,NINCX,IPOLGN

        RETURN

      ELSE
C
C  PARTICLE OUTSIDE STANDARD MESH, IN ADDITIONAL CELL NACELL
C  NRCELL=0
C
        IF (ITFRST == 0) THEN
          ITFRST = 1
          NTETSUR = COUNT(NTBAR(1:4,1:NTET) == 0)
          ALLOCATE (IDROB(NTETSUR,2))
          IOB = 0
          DO ITET = 1,NTET
            DO IS = 1,4
              IF (NTBAR(IS,ITET) == 0) THEN
                IOB = IOB + 1
                IDROB(IOB,:) = (/ ITET,IS /)
              END IF
            END DO
          END DO

          ALLOCATE (ISEEOB(0:NLIMI,NTETSUR/NBITS+1))
          ISEEOB = 0
          ISEEOB(0,:) = -1
          DO I=1,NTETSUR
            IT=IDROB(I,1)
            IS=IDROB(I,2)
            SIG1=SIGN(1,IRICH(1,IS))
            I1=ABS(IRICH(1,IS))
            SIG2=SIGN(1,IRICH(2,IS))
            I2=ABS(IRICH(2,IS))
            A1 = VTETX(I1,IT) * SIG1
            A2 = VTETY(I1,IT) * SIG1
            A3 = VTETZ(I1,IT) * SIG1
            B1 = VTETX(I2,IT) * SIG2
            B2 = VTETY(I2,IT) * SIG2
            B3 = VTETZ(I2,IT) * SIG2
            V1 = A2*B3 - A3*B2
            V2 = A3*B1 - A1*B3
            V3 = A1*B2 - A2*B1
            DO IL=1,NLIMI
              IF (RLB(IL) == 3.) THEN
                S1=(P1(1,IL)-XTETRA(NTECK(ITSIDE(1,IS),IT)))*V1 +
     .             (P1(2,IL)-YTETRA(NTECK(ITSIDE(1,IS),IT)))*V2 +
     .             (P1(3,IL)-ZTETRA(NTECK(ITSIDE(1,IS),IT)))*V3
                S2=(P2(1,IL)-XTETRA(NTECK(ITSIDE(1,IS),IT)))*V1 +
     .             (P2(2,IL)-YTETRA(NTECK(ITSIDE(1,IS),IT)))*V2 +
     .             (P2(3,IL)-ZTETRA(NTECK(ITSIDE(1,IS),IT)))*V3
                S3=(P3(1,IL)-XTETRA(NTECK(ITSIDE(1,IS),IT)))*V1 +
     .             (P3(2,IL)-YTETRA(NTECK(ITSIDE(1,IS),IT)))*V2 +
     .             (P3(3,IL)-ZTETRA(NTECK(ITSIDE(1,IS),IT)))*V3
                IF (ANY( (/ S1,S2,S3 /) > 0.D0)) THEN
                  CALL BITSET (ISEEOB,0,NLIMI,IL,I,1,NBITS)
                ELSE
                  IF (NOPTIM >= NTET) THEN
                    IF (NLIMPB >= NLIMPS) THEN
                      IGJUM3(ITET,IL) = 0
                    ELSE
                      CALL BITSET (IGJUM3,0,NOPTIM,ITET,IL,0,NBITS)
                    END IF
                  END IF
                END IF
              ELSE
                CALL BITSET (ISEEOB,0,NLIMI,IL,I,1,NBITS)
              END IF
            END DO
          END DO
        END IF

        PT=1.D30
        MRSURF=0

        TIMT=1.D30
        NTMZ=0
        NTMS=0
        DO IT=1,NTETSUR
          IF ((MRSURF == 0) .AND.
     .        .NOT.BITGET(ISEEOB,0,NLIMI,MSURF,IT,NBITS)) CYCLE
          I = IDROB(IT,1)
          J = IDROB(IT,2)
          IF (I == IZELLO)  CYCLE
          SIG1=SIGN(1,IRICH(1,J))
          I1=ABS(IRICH(1,J))
          SIG2=SIGN(1,IRICH(2,J))
          I2=ABS(IRICH(2,J))
          A(1:3,1) = (/ VTETX(I1,I), VTETY(I1,I),
     .                  VTETZ(I1,I)/) * SIG1
          A(1:3,2) = (/ VTETX(I2,I), VTETY(I2,I),
     .                  VTETZ(I2,I)/) * SIG2
          A(1:3,3) = (/ -VELX, -VELY, -VELZ /)
          B(1:3) = (/ X0-XTETRA(NTECK(ITSIDE(1,J),I)),
     .                Y0-YTETRA(NTECK(ITSIDE(1,J),I)),
     .                Z0-ZTETRA(NTECK(ITSIDE(1,J),I)) /)
          DET = SARRUS(A)
          IF (ABS(DET) < 1.D-5) CYCLE
          AB = A
          AB(:,1) = B
          XMU = SARRUS(AB)/DET
          AB = A
          AB(:,2) = B
          XETA = SARRUS(AB)/DET
          AB = A
          AB(:,3) = B
          T = SARRUS(AB)/DET
          IF ((XMU .GE.0.D0 .AND. XMU .LE.1.D0) .AND.
     .        (XETA.GE.0.D0 .AND. XETA.LE.1.D0) .AND.
     .        (XMU+XETA <= 1.D0) .AND. (T .GT. 0.)) THEN
            IF (T < TIMT) THEN
              TIMT = T
              NTMZ = I
              NTMS = J
            END IF
          END IF
        END DO
        PT = TIMT
        IZELLO = 0
        IPOLGOO = 0
        IF (NTMZ .NE. 0) THEN
C  SEITENNUMMER DES NEUEN DREIECKS
          IPOLGO = NTMS
          IF (NLTRC) THEN
            WRITE(6,*) 'OUTSIDE; GEHE IN ZELLE ',NTMZ,' TM= ',TM
            WRITE(6,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
          ENDIF
        END IF

        NLSRFX=.FALSE.
        NJUMP=1
        MRSURF=NTMZ
        IPOLGN=NTMS
        NINCX=NTMZ-NRCELL
        IF (NLTRC) WRITE (6,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
     .                           NRCELL,MRSURF,PT,NINCX,IPOLGN
      END IF
      RETURN

12000 CONTINUE
C
C  GENERAL GEOMETRY OPTION: PROVIDE FLIGHT TIME IN CURRENT CELL
C
      IF (ICALL == 0) THEN
        ICALL = 1
        MXSF = MAXVAL(INUMP(1:NSTSI,1))
        IF (MXSF < NSURF) THEN
          NRMSRF = MXSF+1
        ELSE
          DO IRS=1,NSURF
            IF (ALL(INUMP(1:NSTSI,1).NE.IRS)) EXIT
          END DO
          NRMSRF = IRS
          IF (NRMSRF > NSURF) THEN
            WRITE (6,*) ' ERROR IN TIMER, NRMSRF WRONG '
            CALL EXIT
          END IF
        END IF
      END IF
      CALL TIMUSR(NRCELL,X0,Y0,Z0,VELX,VELY,VELZ,NJUMP,
     .            NEWCEL,TIM,ICOS,IERR,NPANU,NLSRFX)
      IF (IERR.NE.0) GOTO 9999
      PT=TIM
      IF (NEWCEL.GT.0) THEN
        NINCX=NEWCEL-NRCELL
        MRSURF=NRMSRF
CPBDR   INMP1I(MRSURF,1,1)=0
      ELSEIF (NEWCEL.LT.0) THEN
        NINCX=ICOS
        MRSURF=INUMP(-NEWCEL,1)
CPBDR   INMP1I(MRSURF,1,1)=-NEWCEL
      ELSEIF (NEWCEL.EQ.0) THEN
        WRITE (6,*) 'NEWCEL=0, EXIT FROM MESH '
        CALL EXIT
      ENDIF
      NLSRFX=.FALSE.
      RETURN

C
9999  CONTINUE
      WRITE (6,*) 'ERROR IN TIMER, EXIT CALLED AT 9999'
      WRITE (6,*) 'NPANU ', NPANU
      CALL EXIT
      END
c === ROUTINE: timet
C
C
      SUBROUTINE TIMET (ZRAD)
C  INPUT
C
C   NRCELL:
C   NTCELL:
C   ZRAD  : DISTANCE (CM) TO THE NEXT RADIAL SURFACE OF 1D STANDARD MESH
C           OR TO NEXT ADDITIONAL SURFACE, TRAVELED IN RADIAL CELL NRCELL
C   PHI   :
C   X01   :
C   Z01   :
C
C  OUTPUT
C
C  IF NLTRZ AND NLTOR
C     NINCZ    :   DIRECTION IN Z-GRID
C     Z01      :
C     NTCELL   :
C     MTSURF   :
C  ELSEIF NLTRA
C    EITHER
C     ISRFCL<3  :   NO ROTATION , NNTCLL=0 , NO PARAMETERS CHANGED
C    OR
C     ISRFCL=3 :   ROTATION CLOCKWISE,
C                  STOP AND RESTART LATER AT MTSURF = NNTCLL
C     ISRFCL=3 :   ROTATION COUNTER CLOCKWISE
C                  STOP AND RESTART LATER AT MTSURF = NNTCLL+1
C
C     PHI      :
C     X01      :
C     Z01      :
C     ZRAD     :   REDUCED TO DISTANCE TO NEXT TOROIDAL SURFACE
C     NNTCLL   :   NEXT POSSIBLE CELL NUMBER IN TOROIDAL MESH
C                  IF PARTICLE TRAVELS REDUCED DISTANCE ZRAD
C     MTSURF   :   NEXT TOROIDAL SURFACE, IF ANY
C     NINCZ    :   DIRECTION IN Z-GRID
C     NINCX    :   RESET TO ZERO
C     MRSURF   :   RESET TO ZERO
C  ELSEIF NLTRP
C    TO BE WRITTEN
C
C  FOR THE TIME BEING: IF NLPOL, REDUCE PATH TO ONE TOROIDAL CELL
C
C
C     CALCULATE TIME SEGMENTS FOR 2D-PROFILES,
C     RADIALLY AND TOROIDALLY RESOLVED
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CGRID
      USE COMPRT
      USE COMSOU
      USE CLGIN

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: ZRAD
      REAL(DP) :: TTT, PHI0, X001, TO, AA, SUM, TU, XTO, EPSTST, XTU,
c slmod begin
c...  Temp, until DUM is removed from PARMMOD.f:
     .          F, FZ, ZRADS, BB, ZRD, Z001, X0TEST, Z0TEST, DZ,
c
c     .          F, FZ, ZRADS, BB, ZRD, DUM, Z001, X0TEST, Z0TEST, DZ,
c slmod end
     .          Y0TEST
      INTEGER :: ITT, ITEST, J2, NERR, IN, J, ICOU, MTTEST, ISTS,
     .           NZSAVE, INCZ, J1, IRSAVE, MTSAVE, LEARCA

      ZRADS=ZRAD
C
      IF (NLTRC) THEN
        CALL LEER(1)
        WRITE (6,*) 'TIMET: ZRAD,NRCELL,NTCELL,MTSURF '
        WRITE (6,*) '      ',ZRAD,NRCELL,NTCELL,MTSURF
        WRITE (6,*) 'INITIAL: X0,X01,Z01,PHI ',X0,X01,Z01,PHI/DEGRAD
      ENDIF
C
      IF (.NOT.NLTRZ) GOTO 1000
C
C  CYLINDRICAL OR CARTHESIAN CO-ORDINATE SYSTEM
C
C  PARTICLE OUTSIDE STANDARD MESH ?
C
      IF (NRCELL.EQ.0) THEN
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT TOROIDAL SURFACES
C
        NCOUT=0
        ZRD=ZRAD
        BB=VELZ+EPS60
        NZSAVE=1
        IF (VELZ.LT.0.D0) NZSAVE=-1
        DO 2903 ISTS=1,NSTSI
          MTTEST=INUMP(ISTS,3)
          IF (MTTEST.NE.0) THEN
C  TEST TOROIDAL SURFACE NO. MTTEST FOR REENTRY
C  TIME FROM Z01 TO ZSURF
            DZ=ZSURF(MTTEST)-Z01
            F=DZ/BB
            IF (NLTRC) WRITE (6,*) 'MTTEST,F,DZ ',MTTEST,F,DZ
            IF (F.LE.ZRD.AND.F.GT.0.D0) THEN
              X0TEST=X0+VELX*F
              Y0TEST=Y0+VELY*F
C  IS THIS RE-ENTRY INSIDE THE STANDARD GRID REGION?
              IF (LEVGEO.EQ.1) THEN
                IF (X0TEST.GE.RSURF(1).AND.X0TEST.LE.RSURF(NR1ST)) THEN
                  IRSAVE=LEARCA(X0TEST,RSURF,1,NR1ST,1,'TIMET 1    ')
                  NCOUT=1
                  JUPC(NCOUT)=1
                  MTSAVE=MTTEST
                  ZRD=F
                ENDIF
              ELSEIF (LEVGEO.GT.1) THEN
                WRITE (6,*) 'ERROR FROM TIMET: '
                WRITE (6,*) 'RE-ENTRY THROUGH TOROIDAL SURFACE NOT '
                WRITE (6,*) 'READY FOR LEVGEO.GT.1 '
                CALL EXIT
              ENDIF
            ENDIF
          ENDIF
2903    CONTINUE
        IF (NCOUT.GT.0.AND.MTSAVE.NE.MTSURF) THEN
C  REENTRY FOUND, REDUCE ZRAD TO F
C  NCOUT=1 AT THIS POINT
          NCOUT=1
          MTSURF=MTSAVE
          NINCZ=NZSAVE
          IRCELL=IRSAVE
          ZRAD=ZRD
          ISRFCL=0
          BLPD(NCOUT)=ZRAD
          MRSURF=0
          MPSURF=0
          MASURF=0
          NINCX=0
          NINCY=0
          IF (NLTRC) THEN
            WRITE (6,*) 'REENTRY FOUND, MTSURF,ZRAD = ',MTSURF,ZRAD
            WRITE (6,*) 'IRCELL ',IRCELL
          ENDIF
          Z01=Z01+VELZ*ZRD
          NTCELL=KUPC(NCOUT)
          ITCELL=NTCELL
          GOTO 5000
        ELSE
C  NO REENTRY FOUND
          NCOUT=1
          KUPC(1)=1
          BLPD(1)=ZRAD
          IF (MRSURF.GT.0) THEN
C  CHECK VALID RANGE ON MRSURF
            Z0TEST=Z00+VELZ*ZRAD
            IF (NLTRC) WRITE (6,*) 'CHECK VALID RANGE: Z0TEST ',Z0TEST
            IF (Z0TEST.GE.ZSURF(1).AND.Z0TEST.LE.ZSURF(NT3RD)) THEN
              ITCELL=LEARCA(Z0TEST,ZSURF,1,NT3RD,1,'TIMET 2     ')
            ELSE
              MRSURF=0
              NINCX=0
            ENDIF
          ENDIF
          MTSURF=0
          NINCZ=0
          IF (NLTRC) THEN
            WRITE (6,*) 'NO REENTRY INTO TOROIDAL GRID FOUND '
          ENDIF
          Z01=Z01+ZRD*VELZ
          GOTO 5000
        ENDIF
C
      ENDIF
C
C  PARTICLE IN STANDARD MESH, RADIAL CELL NRCELL
C
2900  CONTINUE
      Z001=Z01+VELZ*ZRAD
C
      DUM=0.
      NCOUT=1
C
C  J1: CELL INDEX
C  J2  SURFACE INDEX
      J1=NTCELL
      IF (VELZ.LT.0.) THEN
        INCZ=0
        NINCZ=-1
        IF (NLSRFZ) J1=MTSURF-1
      ELSE
        INCZ=1
        NINCZ=1
        IF (NLSRFZ) J1=MTSURF
      ENDIF
      J2=J1+INCZ
C
      NLSRFZ=.FALSE.
C
10    CONTINUE
      IF (J2.LE.0.OR.J2.GT.NT3RD) GOTO 990
C  TIME FROM Z01 TO ZSURF
      DZ=(ZSURF(J2)-Z01)
      BB=VELZ+EPS60
      F=DZ/BB
      IF (F.LE.ZRAD) THEN
        KUPC(NCOUT)=J1
        BLPD(NCOUT)=F-DUM
        DUM=F
C  STOP HISTORY AT NON DEFAULT STANDARD SURFACE J2
        ITEST=INMP3I(0,0,J2)
        IN=ITEST+NLIM
!pb        IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
        IF ((ITEST.NE.0.AND.ILIIN(IN).NE.0).or.(ityp==3)) THEN
          ZRAD=F
          ISRFCL=0
          NTCELL=J1
          MTSURF=J2
          MRSURF=0
          MASURF=0
          NINCX=0
          IPOLGN=0
          ITCELL=NTCELL
          Z01=ZSURF(J2)
          IF (NLTRC) WRITE (6,*) 'STOP AT MTSURF ',MTSURF
          GOTO 5000
        ENDIF
C  NEXT CELL
        J1=J1+NINCZ
        J2=J1+INCZ
        NCOUT=NCOUT+1
        GOTO 10
      ENDIF
C
C  LAST CELL
C
100   CONTINUE
      NTCELL=J1
      MTSURF=0
      NINCZ=0
      KUPC(NCOUT)=J1
      BLPD(NCOUT)=ZRAD-DUM
      ITCELL=NTCELL
      Z01=Z001
C
      GOTO 5000
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN TIMET, Z SURFACE INDEX OUT OF RANGE  '
      WRITE (6,*) 'NPANU,Z0,Z01,ZRAD,VELZ,NTCELL '
      WRITE (6,*)  NPANU,Z0,Z01,ZRAD,VELZ,NTCELL
      WEIGHT=0.
      ZRAD=-1.
      GOTO 5000
C
C  CYLINDER APPROXIMATION FINISHED
C
1000  CONTINUE
C
      IF (.NOT.NLTRA) GOTO 5000
C
C  DISCRETE TOROIDAL APPROXIMATION, ONLY ONE STEP AT A TIME
C
      NERR=0
      NCOUT=1
      KUPC(1)=NTCELL
C
C     IF (NLSRFZ) THEN ....
      NLSRFZ=.FALSE.
c slmod begin (sl)
c...  This may not be required any more.  The problem this addition addressed
c     was related to always launching the target neutral at the segment boundary,
c     which is not done any more.
      NLSRFT=.FALSE.
c slmod end
C
1010  CONTINUE
C  PHI0 IS THE PHI AT THE CENTER OF THE CURRENT TOROIDAL CELL
      PHI0=PHI-ATAN2(Z01,X01)
C
      IF (NLTRC) THEN
        TTT=Z01/(X01*TANAL)
        IF (ABS(TTT).GT.1.+EPS10) THEN
          WRITE (6,*) 'NPANU ',NPANU
          WRITE (6,*) 'X01,Z01 OUT OF RANGE IN TIMET'
          WRITE (6,*) X01,Z01,TTT
          CALL EXIT
        ENDIF
      ENDIF
C
      Z001=Z01+ZRAD*VELZ
      X001=X01+ZRAD*VELX
      IF (ZRAD.LT.1.D30.AND.X01*X001.GT.0.D0) THEN
        ITT=IDINT(REAL(Z001/(X001*TANAL),KIND(1.D0)))
        IF (NLTRC) WRITE (6,*) 'TIMET 1 ',X01,Z01,X001,Z001,ITT
      ELSE
        TO=(Z01-X01*TANAL)/(TANAL*VELX-VELZ)
        XTO=X01+TO*VELX
        TU=(Z01+X01*TANAL)/(-TANAL*VELX-VELZ)
        XTU=X01+TU*VELX
        IF (NLTRC) WRITE (6,*) 'TU,TO ',TU,TO,XTU,XTO
        EPSTST=EPS10*TANAL
        IF (XTO.GT.0..AND.TO.GT.EPSTST) THEN
          ITT=1
        ELSEIF (XTU.GT.0..AND.TU.GT.EPSTST) THEN
          ITT=-1
        ELSE
          ITT=0
          Z001=Z01
          X001=X01
        ENDIF
        IF (NLTRC) WRITE (6,*) 'TIMET 2 ',X01,Z01,ITT
      ENDIF
C
c slmod begin (sl)
      NINCS=0
c slmod end
      IF (ITT) 1100,1200,1300
C
C  NO INTERSECTION WITH TOROIDAL SURFACE
C
1200  CONTINUE
      MTSURF=0
      NNTCLL=IPERID
      Z01=Z001
      X01=X001
C  IN CASE ZRAD=1.D30, IS NEXT STATEMENT IS NONSENSE, BUT CORRECTED
C                      FOR IN SUBR. STDCOL, WHICH MUST BE CALLED NEXT
C                      FOR A POLOIDAL SURFACE (OTHERWISE: ERROR EXIT)
      PHI=PHI0+ATAN2(Z01,X01)
      BLPD(1)=ZRAD
      IF (NLTRC) WRITE (6,*) 'FINAL 1: X01,Z01,PHI ',X01,Z01,PHI/DEGRAD
      GOTO 5000
C
C  PARTICLE LEAVES CELL IN POSITIVE DIRECTION, REDUCE ZRAD
C
1300  CONTINUE
C
C  TIME TO REACH CELL SURFACE: F
C
      AA=Z01-X01*TANAL
      BB=(TANAL*VELX-VELZ)+EPS60
      F=AA/BB
      IF (F.LE.0.) THEN
C  PARTICLE ACCIDENTALLY ON A TOROIDAL SURFACE?
        IF (ABS(AA).LE.EPS10.AND.NERR.LE.1) THEN
          IF (NLTRC) WRITE (6,*) 'TRY AGAIN IN TIMET'
          Z01=Z01-VELZ*EPS10
          X01=X01-VELX*EPS10
          NERR=NERR+1
          GOTO 1010
        ENDIF
        ZRAD=1.D30
        GOTO 9998
      ENDIF
      X01=X01+F*VELX
      Z01=TANAL*X01
      PHI=PHI0+ZHALF
      ZRAD=F
C
      ISRFCL=3
      BLPD(1)=ZRAD
C
      NINCX=0
      NINCZ=1
      IPOLGN=0
      MRSURF=0
      MASURF=0
      NNTCLL=NTCELL+1
c slmod begin (sl)
      NINCS=1
c slmod end
      IF (.NOT.NLTOR) NNTCLL=IPERID+1
      IF (NNTCLL.GE.NTTRA) NNTCLL=1
      MTSURF=NNTCLL
      IF (NLTRC) THEN
        WRITE (6,*) 'FINAL 2: X01,Z01,PHI ',X01,Z01,PHI/DEGRAD
        WRITE (6,*) 'ZRAD,ISRFCL,IPERID ',ZRAD,ISRFCL,IPERID
      ENDIF
      GOTO 5000
C
C  PARTICLE LEAVES CELL IN NEGATIVE DIRECTION, REDUCE ZRAD
C
1100  CONTINUE
C
C  TIME TO REACH CELL SURFACE
C
      AA=Z01+X01*TANAL
      BB=-TANAL*VELX-VELZ+EPS60
      F=AA/BB
      IF (F.LE.0.) THEN
C  PARTICLE ACCIDENTALLY ON A TOROIDAL SURFACE?
        IF (ABS(AA).LE.EPS10.AND.NERR.LE.1) THEN
          IF (NLTRC) WRITE (6,*) 'TRY AGAIN IN TIMET'
          Z01=Z01-VELZ*EPS10
          X01=X01-VELX*EPS10
          NERR=NERR+1
          GOTO 1010
        ENDIF
        ZRAD=1.D30
        GOTO 9998
      ENDIF
      X01=X01+F*VELX
      Z01=-TANAL*X01
      PHI=PHI0-ZHALF
      ZRAD=F
C
      ISRFCL=3
      BLPD(1)=ZRAD
C
      NINCX=0
      NINCZ=-1
      IPOLGN=0
      MRSURF=0
      MASURF=0
      NNTCLL=NTCELL-1
c slmod begin (sl)
c...  NINCS indicates the change in NTRSEG as a result of a 
c     toroidal segment boundary surface:
      NINCS=-1
c slmod end
      IF (.NOT.NLTOR) NNTCLL=IPERID-1
      IF (NNTCLL.LE.0) NNTCLL=NTTRAM
      MTSURF=NNTCLL+1
      IF (NLTRC) THEN
        WRITE (6,*) 'FINAL 3: X01,Z01,PHI ',X01,Z01,PHI/DEGRAD
        WRITE (6,*) 'ZRAD,ISRFCL,IPERID ',ZRAD,ISRFCL,IPERID
      ENDIF
      GOTO 5000
C
C  DISCRETE TOROIDAL APPROXIMATION FINISHED
C
C
5000  CONTINUE
      IF (NLTRC) WRITE (6,*) 'NCOUT= ',NCOUT
      DO 5100 J=1,NCOUT
        CLPD(J)=BLPD(J)
        NUPC(J)=(KUPC(J)-1)*NP2T3
        NCOUNT(J)=KUPC(J)
        IF (CLPD(J).LE.0..OR.KUPC(J).LE.0.OR.
     .      (KUPC(J).GE.NT3RD.AND.NLTOR)) THEN
          WRITE (6,*) 'ERROR DETECTED IN TIMET '
          WRITE (6,*) 'NPANU,J,BLPD,KUPC ',NPANU,J,BLPD(J),KUPC(J)
        ENDIF
        IF (NLTRC) THEN
          WRITE (6,*) 'TIMET: J,BLPD,NUPC,NCOUNT ',
     .                        J,BLPD(J),NUPC(J),NCOUNT(J)
        ENDIF
5100  CONTINUE
      IF (NLTRC) THEN
        WRITE (6,*) 'MTSURF,NLSRFZ,NINCZ,IRCELL,IPCELL '
        WRITE (6,*)  MTSURF,NLSRFZ,NINCZ,IRCELL,IPCELL
        IF (NLTOR) WRITE (6,*) 'INMP3I ',INMP3I(IRCELL,IPCELL,MTSURF)
      ENDIF
C
      NCOU=NCOUT
C
      SUM=0.
      DO 5110 ICOU=1,NCOU
        SUM=SUM+CLPD(ICOU)
C       WRITE (6,*) 'ICOU,NCOUNT,CLPD ',ICOU,NCOUNT(ICOU),CLPD(ICOU)
5110  CONTINUE
      IF (MTSURF.EQ.0.AND.ABS(SUM-ZRADS).GT.EPS10) THEN
        WRITE (6,*) 'ERROR IN TIMET: NPANU,SUM,ZRADS ',
     .                               NPANU,SUM,ZRADS
        WRITE (6,*) 'TRY TO KILL PARTICLE ASAP '
        SUM=-1.0D0
      ENDIF
      ZRAD=SUM
C
      RETURN
C
9998  CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN TIMET. RETURN ZRAD=1.D30'
      WRITE (6,*) 'NPANU,AA,BB ',NPANU,AA,BB
      RETURN
9999  CONTINUE
      WRITE (6,*) 'INVALID OPTION IN TIMET. EXIT CALLED '
      CALL EXIT
      END
c === ROUTINE: xshadd
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
c === ROUTINE: yshadd
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
c === ROUTINE: zshadd
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
