C@PROCESS IL(DIM) OPT(3) NOSDUMP NOGOSTMT
      SUBROUTINE GRVAR(X,Y,VAR,N,KLIP)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C
C     GRVAR zeichnet N Fehlerbalken an die Stellen X(I),Y(I) von
C     Y(i)-VAR(i) bis Y(i)+VAR(i)           bei KLIP=0,    von
C     X(i)-VAR(i) bis X(i)+VAR(i)           bei KLIP^=0.
C     Die Groesse der Punkte & End-Querstriche wird mit GRMRKS bestimmt.
C     Die Liniendicke mit GRSPTS, Farbe mit GRNWPN. Clipping moeglich.
C
C***********************************************************************
      IMPLICIT NONE
C
      INTEGER N,KLIP
      REAL X(N),Y(N),VAR(N),FAC
      PARAMETER (FAC=.5)

      REAL PP,SHI,SYSIZ,CSIZ,CANG
      INTEGER IFONT,INTLIN,ICOL
      LOGICAL ROT
      COMMON /GRPP/ PP(8),IFONT,SHI(3),INTLIN,CSIZ,CANG,ICOL,SYSIZ,ROT
CDEC$ PSECT /GRPP/ NOSHR
      SAVE /GRPP/

      REAL STX(7),STY(7),B1X(2),B1Y(2)
      REAL B2X(2),B2Y(2),FAX,FAY,SYX,SYY,SYXF,SYYF
      INTEGER I

      FAX = (PP(7)-PP(5))/(PP(3)-PP(1))
      FAY = (PP(8)-PP(6))/(PP(4)-PP(2))
      SYX = SYSIZ*FAX*.5
      SYY = SYSIZ*FAY*.5
      SYXF = SYX*FAC
      SYYF = SYY*FAC

      IF (KLIP.EQ.0) THEN
         DO 1 I=1,N
            B1X(1) = X(I)-SYX
            B1Y(1) = Y(I)+VAR(I)
            B1X(2) = X(I)+SYX
            B1Y(2) = B1Y(1)
            B2X(1) = B1X(1)
            B2Y(1) = Y(I)-VAR(I)
            B2X(2) = B1X(2)
            B2Y(2) = B2Y(1)
            STX(1) = X(I)
            STY(1) = B1Y(1)
            STX(2) = X(I)
            STY(2) = Y(I)+SYYF
            STX(3) = X(I)-SYXF
            STY(3) = Y(I)
            STX(4) = X(I)
            STY(4) = Y(I)-SYYF
            STX(5) = X(I)+SYXF
            STY(5) = Y(I)
            STX(6) = X(I)
            STY(6) = STY(2)
            STX(7) = X(I)
            STY(7) = B2Y(1)
            CALL GRLN(B1X,B1Y,2)
            CALL GRLN(B2X,B2Y,2)
            CALL GRLN(STX,STY,7)
    1    CONTINUE
      ELSE
         DO 2 I=1,N
            B1Y(1) = Y(I)-SYY
            B1X(1) = X(I)-VAR(I)
            B1Y(2) = Y(I)+SYY
            B1X(2) = B1X(1)
            B2Y(1) = B1Y(1)
            B2X(1) = X(I)+VAR(I)
            B2Y(2) = B1Y(2)
            B2X(2) = B2X(1)
            STY(1) = Y(I)
            STX(1) = B1X(1)
            STY(2) = Y(I)
            STX(2) = X(I)-SYXF
            STY(3) = Y(I)-SYYF
            STX(3) = X(I)
            STY(4) = Y(I)
            STX(4) = X(I)+SYXF
            STY(5) = Y(I)+SYYF
            STX(5) = X(I)
            STY(6) = Y(I)
            STX(6) = STX(2)
            STY(7) = Y(I)
            STX(7) = B2X(1)
            CALL GRLN(B1X,B1Y,2)
            CALL GRLN(B2X,B2Y,2)
            CALL GRLN(STX,STY,7)
    2    CONTINUE
      ENDIF

      END
