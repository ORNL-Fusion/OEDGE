C@process OPT(3) nogostmt nosdump
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE 12.6.89 GROTEN  STRICHDICKE DURCH PISY+2000
C UPDATE 16.8.90 GROTEN  STRICHDICKE DURCH PISY+2000*I, NEUES GKS, CRAY
C UPDATE 25.9.90 GROTEN
C UPDATE 23.1.92 Busch Vektor SYMVEC mit Symbolzuordnung analog GRBLD
C                eingefuegt
      SUBROUTINE GRLGND(RT)
      LOGICAL   B100,BSTR,W(14),V,SHORT
      CHARACTER(len=97) TX
      CHARACTER(len=96)FTX
      integer pisy, symvec(13)
      REAL   X(14),
     F       Y(14),
     F       SH(3,5)
      CHARACTER(len=8) RT(6,14),FLOAT,F(12)
      EQUIVALENCE  (F(1),FTX)
      COMMON /GRCIL/ CILBER(8)
CDEC$ PSECT /GRCIL/ NOSHR
C---- Grisy set in GRBLD
      COMMON /GRISY/ PISY(14)
CDEC$ PSECT /GRISY/ NOSHR
      SAVE /GRCIL/,/GRISY/

      DATA W/6*.FALSE.,6*.TRUE.,2*.FALSE./
      DATA   X
     F     /3.2,17.6,3.2,17.6,3.2,17.6,-2.5,-1.8,-1.1,3*0.,3.2,17.6/
      DATA   Y
     F     /-1.4,-1.4,-2.2,-2.2,-3.,-3.,8*3.5/
      DATA   SH     / .25   , .25   , .25   ,
     F                .25   , .0625 , .03125,
     F                .03125, .09375, .03125,
     F                .5    , .25   , .25   ,
     F                1.     , .0625 ,1.      /
      DATA SYMVEC/2,3,4,5,-101,-102,-103,-108,-109,-110,-111,-112,-113/

      L=0
      XX=CILBER(1)
      YY=CILBER(2)
      X(10)=XX+1.2
      X(11)=XX+2.
      X(12)=XX+2.8
      Y(13)=YY+.5
      Y(14)=Y(13)
      SHORT=.TRUE.
      DO 13 I=1,14
      IF (SHORT) GOTO 20
      SHORT=.TRUE.
      GOTO 13
   20 XX=X(I)
      YY=Y(I)
      M=0
      IF (.NOT.(I.EQ.1 .OR. I.EQ.3 .OR. I.EQ.5 .OR. I.EQ.13)) GOTO 21
      II=12
      J=96
      SHORT=.FALSE.
      GOTO 22
   21 II=6
      J=48
   22 FLOAT=RT(1,I)
      V=W(I)
      R=0.
      IF (V) R=90.
      IF (FLOAT(1:1).NE.'#') GOTO 10
      M=1
      IF (SHORT) GOTO 23
      J=86
      II=11
      GOTO 24
   23 J=38
      II=5
   24 L=L+1
      K=PISY(L)
C
C FARBE, ZEICHENDICKE, SYMBOLGROESSE
C
      jj=mod(k,2000)
      IFARB=jj/200+1
      IF (K.LT.2000) GOTO 1
      N=18+K/2000*2
      G=.5
      K=MOD(K,2000)
      GOTO 4
    1 N=18
      G=.3
    4 CALL GRSPTS(N)
      CALL GRNWPN(IFARB)
      k=mod(k,200)
C
C  VERLAUF DER STRECKE, START FUER TEXT
C
      IF (V) GOTO 5
      X1=XX+.3
      Y1=YY+.15
      X5=XX+2.5
      Y5=YY+.15
      XX=XX+3.
      GOTO 6
    5 X1=XX-.15
      Y1=YY+.3
      X5=XX-.15
      Y5=YY+2.5
      YY=YY+3.
C
C  STRICHELUNG ?
C
    6 BSTR=K.GT.113
      B100=K.GE.100
      IF (B100) K=K-100
      DX=(X5-X1)*.5
      DY=(Y5-Y1)*.5
      IF (K.LE.13) GOTO 7
      IF (BSTR) CALL GRDSH(SH(1,K-13),SH(2,K-13),SH(3,K-13))
      CALL GRJMP(X1,Y1)
      CALL GRDRW(X5,Y5)
      IF (BSTR) CALL GRDSH(1.,0.,1.)
      GOTO 10
    7 IF(K.GE.0) GOTO 71
      CALL GRJMP(X1,Y1)
      CALL GRDRW(X1+.02,Y1+.02)
      GOTO 72
   71 CALL GRCHRC(G,R,N)
      ISYM=SYMVEC(K)
      CALL GRJMPS(X1,Y1,ISYM)
   72 DO 9 N=1,2
      X1=X1+DX
      Y1=Y1+DY
      IF (B100) GOTO 8
      IF(K.LE.0) GOTO 73
      ISYM=SYMVEC(K)
      CALL GRJMPS(X1,Y1,ISYM)
      GOTO 9
   73 CALL GRJMP(X1,Y1)
      CALL GRDRW(X1+.02,Y1+.02)
      GOTO 9
    8 ISYM=SYMVEC(K)
      CALL GRDRWS(X1,Y1,ISYM)
    9 CONTINUE
10    CALL GRNWPN(1)
      CALL GRCHRC(.3,R,18)
      DO 11 N=1,II
        F(N)=RT(N,I)
 11   CONTINUE
      tx(:j)=ftx(m+1:m+j)
      TX(J+1:j+1)='@'
      CALL GRTXT(XX,YY,J,TX)
   13 CONTINUE
      CALL GRSPTS(18)
      END
