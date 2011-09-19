C
C
      SUBROUTINE MSHADJ (X1,Y1,X2,Y2,XPLG,YPLG,NPLP,NPOINT,M1,NDX,IR)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1(*),Y1(*),X2(*),Y2(*)
      INTEGER, INTENT(IN) ::M1, NDX, IR
      REAL(DP), INTENT(INOUT) :: XPLG(M1,*),YPLG(M1,*)
      INTEGER, INTENT(INOUT) :: NPOINT(2,*), NPLP
      REAL(DP) :: EPS, D, D1
      INTEGER :: NPRT, NP, I
      LOGICAL :: LDUMCL,LPRT
C  LPRT=.TRUE.  : VALID PART
C  LDUMCL=.TRUE.: DUMMY CELL
C
      EPS=1.E-20
      NPRT=0
      NP=0
      LDUMCL=.FALSE.
      LPRT=.FALSE.
C
      DO 1 I=1,NDX
        D=ABS((X1(I)-X2(I))**2+(Y1(I)-Y2(I))**2)
        IF (D.LE.EPS) THEN
          IF (LPRT) THEN
C  ENDE DER VALID ZELLEN
            NPOINT(2,NPRT)=NP
            LPRT=.FALSE.
            LDUMCL=.TRUE.
          ENDIF
        ELSEIF (D.GT.EPS) THEN
C  STARTPUNKT DIESES POLYGONS = ERSTER VALID PUNKT?
          IF (.NOT.LPRT.AND..NOT.LDUMCL) THEN
            LPRT=.TRUE.
            NPRT=NPRT+1
            NP=NP+1
            XPLG(IR,NP)=X1(I)
            YPLG(IR,NP)=Y1(I)
            NPOINT(1,NPRT)=NP
          ELSEIF (.NOT.LPRT.AND.LDUMCL) THEN
            D1=ABS((X1(I+1)-X2(I+1))**2+(Y1(I+1)-Y2(I+1))**2)
            IF (D1.GT.EPS) THEN
C  ENDE DER DUMMYZELLEN: SUCHE NAECHSTE ZELLE MIT D1.GT.0.
              LDUMCL=.FALSE.
              LPRT=.TRUE.
              NPRT=NPRT+1
              NPOINT(1,NPRT)=NP
            ENDIF
          ENDIF
        ENDIF
        NP=NP+1
        XPLG(IR,NP)=X2(I)
        YPLG(IR,NP)=Y2(I)
1     CONTINUE
      NPOINT(2,NPRT)=NP
      NPLP=NPRT
      RETURN
      END
