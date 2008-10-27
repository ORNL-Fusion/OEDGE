C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS TB04AD
      SUBROUTINE TB04AD(N,X,F,D,A)
C STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
      DOUBLE PRECISION A,D,F,H1,H2,P,X
      DIMENSION X(N),F(N),D(N),A(*)
      DATA NP/6/
C F(I) ARE THE FUNCTION VALUES AT THE POINTS X(I) FOR I=1,N AND
C THE SPLINE DERIVATIVES D(I) ARE FOUND.  THE DIMENSION OF A MUST
C NOT BE LESS THAN 3*N. PERIPHERAL NP MUST BE AN OUTPUT MEDIUM.
      DO 5 I=2,N
      IF(X(I)-X(I-1))1,1,5
1     WRITE(NP,3)I
3      FORMAT(29H RETURN FROM TB04AD BECAUSE X,I3,13H OUT OF ORDER)
      A(1)=1.0D0
      RETURN
5     CONTINUE
      DO 30 I=1,N
      J=2
      IF(I-1)6,10,6
6     J=N-1
      IF(I.EQ.N)GO TO 10
      H1=1.0D0/(X(I)-X(I-1))
      H2=1.0D0/(X(I+1)-X(I))
      A(3*I-2)=H1
      A(3*I-1)=2.0D0*(H1+H2)
      A(3*I)=H2
      D(I)=3.0D0*(F(I+1)*H2*H2+F(I)*(H1*H1-H2*H2)-F(I-1)*H1*H1)
      GO TO 30
10    H1=1.0D0/(X(J)-X(J-1))
      H2=1.0D0/(X(J+1)-X(J))
      A(3*I-2)=H1*H1
      A(3*I-1)=H1*H1-H2*H2
      A(3*I)=-H2*H2
      D(I)=2.D0*(H1*H1*H1*(F(J)-F(J-1))+H2*H2*H2*(F(J)-F(J+1)))
30    CONTINUE
      P=A(4)/A(1)
      A(5)=A(5)-P*A(2)
      A(6)=A(6)-P*A(3)
      D(2)=D(2)-P*D(1)
      DO 50 I=3,N
      K=3*I-4
      P=A(K+2)/A(K)
      A(K+3)=A(K+3)-P*A(K+1)
      D(I)=D(I)-P*D(I-1)
      IF(I.NE.N-1)GO TO 50
      P=A(K+5)/A(K)
      A(K+5)=A(K+6)-P*A(K+1)
      A(K+6)=A(K+7)
      D(N)=D(N)-P*D(N-2)
50    CONTINUE
      D(N)=D(N)/A(3*N-1)
      DO 60 I=3,N
      J=N+2-I
60    D(J)=(D(J)-A(3*J)*D(J+1))/A(3*J-1)
      D(1)=(D(1)-D(2)*A(2)-D(3)*A(3))/A(1)
      A(1)=0.0D0
      RETURN
      END
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS TG01BD
      DOUBLE PRECISION FUNCTION TG01BD(IX,N,U,S,D,X)
C
C**********************************************************************
C
C      TG01BD - FUNCTION ROUTINE TO EVALUATE A CUBIC SPLINE GIVEN SPLINE
C     VALUES AND FIRST DERIVATIVE VALUES AT THE GIVEN KNOTS.
C
C     THE SPLINE VALUE IS DEFINED AS ZERO OUTSIDE THE KNOT RANGE,WHICH
C     IS EXTENDED BY A ROUNDING ERROR FOR THE PURPOSE.
C
C                 F = TG01BD(IX,N,U,S,D,X)
C
C       IX    ALLOWS CALLER TO TAKE ADVANTAGE OF SPLINE PARAMETERS SET
C             ON A PREVIOUS CALL IN CASES WHEN X POINT FOLLOWS PREVIOUS
C             X POINT. IF IX < 0 THE WHOLE RANGE IS SEARCHED FOR KNOT
C             INTERVAL; IF IX > 0 IT IS ASSUMED THAT X IS GREATER THAN
C             THE X OF THE PREVIOUS CALL AND SEARCH STARTED FROM THERE.
C       N     NUMBER OF KNOTS.
C       U     THE KNOTS.
C       S     THE SPLINE VALUES.
C       D     THE FIRST DERIVATIVE VALUES OF THE SPLINE AT THE KNOTS.
C       X     THE POINT AT WHICH THE SPLINE VALUE IS REQUIRED.
C       F     THE VALUE OF THE SPLINE AT THE POINT X.
C
C                                      MODIFIED JULY 1970
C
C**********************************************************************
C
      DOUBLE PRECISION A,B,D,H,Q1,Q2,S,SS,U,X,Z
C
C     ALLOWABLE ROUNDING ERROR ON POINTS AT EXTREAMS OF KNOT RANGE
C     IS 2**IEPS*MAX(!U(1)!,!U(N)!).
      DIMENSION U(1),S(1),D(1)
      DATA IFLG/0/,IEPS/-50/
C
C       TEST WETHER POINT IN RANGE.
      IF(X.LT.U(1)) GO TO 990
      IF(X.GT.U(N)) GO TO 991
C
C       JUMP IF KNOT INTERVAL REQUIRES RANDOM SEARCH.
      IF(IX.LT.0.OR.IFLG.EQ.0) GO TO 12
C       JUMP IF KNOT INTERVAL SAME AS LAST TIME.
      IF(X.LE.U(J+1)) GO TO 8
C       LOOP TILL INTERVAL FOUND.
    1 J=J+1
   11 IF(X.GT.U(J+1)) GO TO 1
      GO TO 7
C
C       ESTIMATE KNOT INTERVAL BY ASSUMING EQUALLY SPACED KNOTS.
   12 J=DABS(X-U(1))/(U(N)-U(1))*(N-1)+1
C       ENSURE CASE X=U(N) GIVES J=N-1.
      J=MIN0(J,N-1)
C       INDICATE THAT KNOT INTERVAL INSIDE RANGE HAS BEEN USED.
      IFLG=1
C       SEARCH FOR KNOT INTERVAL CONTAINING X.
      IF(X.GE.U(J)) GO TO 11
    2 J=J-1
      IF(X.LT.U(J)) GO TO 2
C
C       CALCULATE SPLINE PARAMETERS FOR JTH INTERVAL.
    7 H=U(J+1)-U(J)
      Q1=H*D(J)
      Q2=H*D(J+1)
      SS=S(J+1)-S(J)
      B=3D0*SS-2D0*Q1-Q2
      A=Q1+Q2-2D0*SS
C
C       CALCULATE SPLINE VALUE.
    8 Z=(X-U(J))/H
      TG01BD=((A*Z+B)*Z+Q1)*Z+S(J)
      RETURN
C       TEST IF X WITHIN ROUNDING ERROR OF U(1).
  990 IF(X.LE.U(1)-2D0**IEPS*DMAX1(DABS(U(1)),DABS(U(N)))) GO TO 99
      J=1
      GO TO 7
C       TEST IF X WITHIN ROUNDING ERROR OF U(N).
  991 IF(X.GE.U(N)+2D0**IEPS*DMAX1(DABS(U(1)),DABS(U(N)))) GO TO 99
      J=N-1
      GO TO 7
   99 IFLG=0
C       FUNCTION VALUE SET TO ZERO FOR POINTS OUTSIDE THE RANGE.
      TG01BD=0D0
      RETURN
      END
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS VC03AD
C###### CALLS   VB06
      SUBROUTINE VC03AD (M,N,XD,YD,WD,RD,XN,FN,GN,DN,THETA,IPRINT,W)
      DOUBLE PRECISION DN,FN,GN,HLF,HS,PRP,RD,RP,SA,SR,SW,THETA,
     *                 TMAX,W,WD,XD,XN,YD
      DIMENSION XD(1),YD(1),WD(1),RD(1),XN(*),FN(*),GN(1),DN(1),
     1THETA(1),W(1)
      NDIMAX=N
C     OMIT DATA WITH ZERO WEIGHTS
      MM=1
      J=1
      HLF=0.5D0
      DO 1 I=2,M
      IF (WD(I)) 3,2,3
    2 IF (I-M) 4,3,3
    4 W(J)=XD(I)
      W(J+1)=YD(I)
      W(J+2)=DFLOAT(I)
      J=J+3
      GO TO 1
    3 MM=MM+1
      XD(MM)=XD(I)
      YD(MM)=YD(I)
      WD(MM)=WD(I)
    1 CONTINUE
      J=1
      K=MM
    7 IF (K-M) 5,6,6
    5 K=K+1
      XD(K)=W(J)
      YD(K)=W(J+1)
      WD(K)=W(J+2)
      J=J+3
      GO TO 7
C     INITIALIZATION OF ITERATIONS
    6 JA=1
      IF (WD(1)) 9,8,9
    8 JA=2
    9 N=5
      IF(N-NDIMAX)100,100,110
  100 XN(1)=XD(1)
      XN(5)=XD(MM)
      XN(3)=HLF*(XN(1)+XN(5))
      XN(2)=HLF*(XN(1)+XN(3))
      XN(4)=HLF*(XN(3)+XN(5))
      IW=7*MM+46
      DO 10 I=1,5
      J=IW+I
      W(J)=XN(I)
   10 CONTINUE
      IP=-IPRINT
C     CALCULATE THE HISTOGRAMS FOR NEW SCALE FACTORS
   11 J=0
      K=1
      SA=0.0D0
   12 SA=SA+WD(K)**2
      K=K+1
      IF (XD(K)-XD(1)) 12,12,13
   13 SA=(SA+HLF*WD(K)**2)/(XD(K)-XD(1))
   14 J=J+1
      IF (XD(K)-XN(J+1)) 15,15,16
   16 GN(J)=SA*(XN(J+1)-XN(J))
      GO TO 14
   15 GN(J)=SA*(XD(K)-XN(J))+HLF*WD(K)**2
   17 K=K+1
      IF (K-MM) 18,18,19
   18 IF (XD(K)-XN(J+1)) 20,20,21
   20 GN(J)=GN(J)+WD(K)**2
      GO TO 17
   21 SA=HLF*(WD(K-1)**2+WD(K)**2)/(XD(K)-XD(K-1))
      GN(J)=GN(J)-HLF*WD(K-1)**2+SA*(XN(J+1)-XD(K-1))
      GO TO 14
C     CALCULATE THE NEW SCALE FACTORS
   19 K=IW+2
      GN(1)=0.00025216D0*GN(1)*(W(K)-W(K-1))**8/(XN(2)-XN(1))
      DO 22 J=3,N
      IF (XN(J)-W(K)) 23,23,24
   24 K=K+1
   23 GN(J-1)=0.00025216D0*GN(J-1)*(W(K)-W(K-1))**8/(XN(J)-XN(J-1))
      GN(J-2)=DLOG(GN(J-2)+GN(J-1))
   22 CONTINUE
      NN=N-2
      HS=1.386294D0
      DO 97 J=2,NN
      GN(J)=DMIN1(GN(J),GN(J-1)+HS)
   97 CONTINUE
      J=NN-1
   98 GN(J)=DMIN1(GN(J),GN(J+1)+HS)
      J=J-1
      IF (J) 99,99,98
   99 DO 25 J=3,N
      THETA(J-1)=DSQRT(DEXP(GN(J-2))/(XN(J)-XN(J-2)))
   25 CONTINUE
C     CALCULATE THE SPLINE APPROXIMATION WITH CURRENT KNOTS
      CALL VB06AD (MM,N,XD,YD,WD,RD,XN,FN,GN,DN,THETA,IP,W)
C     APPLY STATISTICAL TEST FOR EXTRA KNOTS
      J=IW+1
      JJ=0
      K=JA
      TMAX=0.0D0
      IIS=1
   26 KC=-1
      SW=0.0D0
      SR=0.0D0
      RP=0.0D0
      J=J+1
      JJ=JJ+1
      W(JJ)=0.0D0
   27 IF (W(J)-XD(K)) 28,29,30
   30 KC=KC+1
      SW=SW+RD(K)**2
      SR=SR+RP*RD(K)
      RP=RD(K)
      K=K+1
      GO TO 27
   29 IF (WD(K)) 31,28,31
   31 KC=KC+1
      SW=SW+RD(K)**2
      SR=SR+RP*RD(K)
   28 IF (SR) 32,32,33
   33 SW=(SW/SR)**2
      RP=DSQRT(SR/DFLOAT(KC))
      IF (DFLOAT(KC)-SW) 32,32,34
   34 GO TO (35,36,37),IIS
   35 PRP=RP
      IIS=3
      IF (DFLOAT(KC)-2.D0*SW) 38,38,39
   37 W(JJ-1)=PRP
      TMAX=DMAX1(TMAX,PRP)
   39 IIS=2
   36 W(JJ)=RP
      TMAX=DMAX1(TMAX,RP)
      GO TO 38
   32 IIS=1
   38 IF (W(J)-XN(N)) 26,40,40
C     TEST WHETHER ANOTHER ITERATION IS REQUIRED
   40 IF (TMAX) 41,41,42
C     CALCULATE NEW TREND ARRAY, INCLUDING LARGER TRENDS ONLY
   42 TMAX=HLF*TMAX
      I=0
      J=1
      JW=1
      K=IW+1
      THETA(JW)=W(K)
   43 I=I+1
      K=K+1
      IF (W(I)-TMAX) 44,44,45
   44 JW=JW+1
      THETA(JW)=W(K)
   46 FN(J)=0.0D0
      J=J+1
      IF (W(K)-XN(J)) 47,47,46
   45 JW=JW+2
      THETA(JW-1)=HLF*(W(K-1)+W(K))
      THETA(JW)=W(K)
      IF (XN(J+1)-THETA(JW-1)) 46,46,48
   48 FN(J)=1.0D0
      J=J+1
   47 IF (J-N) 43,49,49
C     MAKE KNOT SPACINGS BE USED FOUR TIMES
   49 IK=1
      KL=1
      FN(2)=DMAX1(FN(1),FN(2))
      GO TO 102
   50 K=KL+3
   51 IF (FN(K)) 52,52,53
   52 K=K-1
      IF (K-KL) 74,74,51
   53 K=K-1
      FN(K)=1.0D0
      IF (K-KL) 74,74,53
  102 K=KL+3
   54 K=K+1
      IF (K-N) 55,56,56
   55 IF (XN(K+1)-XN(K)-1.5*(XN(K)-XN(K-1))) 54,54,56
   56 KU=K
      FN(K-2)=DMAX1(FN(K-2),FN(K-1))
   57 KKU=K
   58 K=K-1
      IF (K-KL) 59,59,60
   60 IF (XN(K)-XN(K-1)-1.5*(XN(K+1)-XN(K))) 58,58,61
   61 FN(K+1)=DMAX1(FN(K),FN(K+1))
   59 KKL=K
      KZ=4
      IF (FN(K)) 62,62,63
   63 K=K+1
      IF (K-KKU) 64,65,65
   64 IF (FN(K)) 66,66,63
   66 KZ=0
   62 KZ=KZ+1
      K=K+1
      IF (K-KKU) 67,65,65
   67 IF (FN(K)) 62,62,68
   68 IF (KZ-3) 69,69,70
   69 J=K-KZ
   71 FN(J)=1.0D0
      J=J+1
      IF (J-K) 71,63,63
   70 IF (K+1-KKU) 72,65,65
   72 K=K+1
      FN(K)=1.0D0
      GO TO 63
   65 IF (KL-KKL) 73,50,50
   73 FN(KKL-2)=DMAX1(FN(KKL-2),FN(KKL+1))
      FN(KKL-1)=DMAX1(FN(KKL-1),FN(KKL+3))
   75 K=KKL-4
   78 IF (FN(K)) 76,76,77
   76 K=K+1
      IF (K-KKL) 78,79,79
   77 FN(K)=1.0D0
      K=K+1
      IF (K-KKL) 77,79,79
   79 GO TO (57,80),IK
   74 IF (KU-N) 81,82,82
   81 KL=KU
      FN(KL+1)=DMAX1(FN(KL+1),FN(KL-2))
      FN(KL)=DMAX1(FN(KL),FN(KL-4))
      GO TO 102
   82 IK=2
      KKL=N
      GO TO 75
C     INSERT EXTRA KNOTS FOR NEW APPROXIMATION
   80 DO 83 J=1,N
      GN(J)=XN(J)
   83 CONTINUE
      NN=1
      DO 84 J=2,N
      IF (FN(J-1)) 85,85,86
   86 NN=NN+1
      XN(NN)=HLF*(GN(J-1)+GN(J))
   85 NN=NN+1
      XN(NN)=GN(J)
   84 CONTINUE
      IF(N-NDIMAX)101,110,110
  110 PRINT 111,N
  111 FORMAT(///' ARRAY SIZES TOO SMALL.  N =',I6,///)
      N=-N
      GO TO 90
  101 N=NN
      IW=7*MM+8*N+6
      DO 87 J=1,JW
      I=IW+J
      W(I)=THETA(J)
   87 CONTINUE
      GO TO 11
C     RESTORE DATA WITH ZERO WEIGHTS
   41 IF (MM-M) 88,89,89
   89 IF (IPRINT) 90,90,91
   88 J=-2
      K=MM
   92 J=J+3
      K=K+1
      W(J)=XD(K)
      W(J+1)=YD(K)
      W(J+2)=WD(K)
      IF (K-M) 92,93,93
   93 I=W(J+2)+HLF
   94 IF (K-I) 95,95,96
   96 XD(K)=XD(MM)
      YD(K)=YD(MM)
      WD(K)=WD(MM)
      K=K-1
      MM=MM-1
      GO TO 94
   95 XD(K)=W(J)
      YD(K)=W(J+1)
      WD(K)=0.0D0
      K=K-1
      J=J-3
      IF (J) 91,91,93
   91 CALL VB06AD (M,N,XD,YD,WD,RD,XN,FN,GN,DN,THETA,IABS(IPRINT),W)
   90 RETURN
      END
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS VB06AD
      SUBROUTINE VB06AD(M,N,XD,YD,WD,RD,XN,FN,GN,DN,THETA,IPRINT,W)
      DOUBLE PRECISION ALPHA,BETA,C,DN,DSTAR,FN,GAMMA,GN,H,PREC
      DOUBLE PRECISION RD,THETA,W,WD,XD,XN,YD,ZERO,ONE,TWO,SIX
      LOGICAL SWITCH
      DIMENSION XD(1),YD(1),WD(1),RD(1),XN(1),FN(1),GN(1),DN(1),
     1THETA(1),W(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,SIX/6.0D0/
C     SET THE RELATIVE ACCURACY OF THE COMPUTER ARITHMETIC
      PREC=16.0D0**(-14)
C     RESERVE THE FIRST FIFTEEN LOCATIONS OF W FOR
C     COMPUTED VALUES AND THIRD DERIVATIVES OF B-SPLINES
      W(1)=ZERO
      W(7)=ZERO
C     THEN SET THE KNOT POSITIONS INCLUDING SIX OUTSIDE THE RANGE
C     AND DELETE ANY KNOTS NOT IN STRICTLY ASCENDING ORDER
      I=1
      NN=1
   10 W(NN+18)=XN(I)
   20 IF (I.GE.N) GO TO 30
      I=I+1
      IF (XN(I).LE.W(NN+18)) GO TO 20
      NN=NN+1
      XN(NN)=XN(I)
      THETA(NN)=THETA(I)
      GO TO 10
   30 IF (NN.GE.N) GO TO 50
      N=NN
      PRINT 40
   40 FORMAT (//5X,'SOME KNOTS HAVE BEEN DELETED BY VB06AD')
   50 IF (NN.GE.2) GO TO 70
      PRINT 60
   60 FORMAT (//5X,'THE EXECUTION OF VB06AD HAS BEEN STOPPED ',
     1'BECAUSE THERE ARE TOO FEW KNOTS')
      GO TO 770
   70 DO 80 I=1,3
      W(19-I)=TWO*W(20-I)-W(21-I)
   80 W(NN+I+18)=TWO*W(NN+I+17)-W(NN+I+16)
C     DELETE ANY DATA POINTS THAT ARE NOT IN ASCENDING ORDER
      I=1
      MM=1
      XD(MM)=DMIN1(DMAX1(XD(1),XN(1)),XN(NN))
   90 IF (I.GE.M) GO TO 100
      I=I+1
      IF (XD(I).LT.XD(MM)) GO TO 90
      MM=MM+1
      XD(MM)=DMIN1(XD(I),XN(NN))
      YD(MM)=YD(I)
      WD(MM)=WD(I)
      GO TO 90
  100 IF (MM.GE.M) GO TO 120
      M=MM
      PRINT 110
  110 FORMAT (//5X,'SOME DATA POINTS HAVE BEEN DELETED BY VB06AD')
C     COMPUTE THE THIRD DERIVATIVES OF THE B-SPLINES AT XN(1)
  120 SWITCH=.TRUE.
      IN=1
      GO TO 280
  130 DO 140 I=2,5
      W(I+6)=W(I)
  140 W(I+10)=W(I)
C     INITIALIZE THE FACTORIZATION OF THE LEAST SQUARES MATRIX
      KS=NN+22
      KU=KS+22
      DO 150 I=KS,KU
  150 W(I)=ZERO
      IM=1
      IN=2
      KK=KS
      SWITCH=.FALSE.
C     TEST WHETHER A DATA POINT OR A KNOT DEFINES THE NEXT EQUATION
  160 KR=KK
      NC=5
      IF (IM.GT.MM) GO TO 260
      IF (XD(IM).GT.XN(IN)) GO TO 260
C     CALCULATE THE VALUES OF THE B-SPLINES AT XD(IM)
      W(2)=ONE/(W(IN+18)-W(IN+17))
      DO 190 J=2,4
      I=J+1
      W(I)=ZERO
  180 II=IN+I+16
      C=(XD(IM)-W(II-J))*W(I-1)+(W(II)-XD(IM))*W(I)
      IF (DABS(C).LE.PREC) C=0.
      W(I)=C/(W(II)-W(II-J))
      I=I-1
      IF (I.GE.2) GO TO 180
  190 CONTINUE
C     COMPLETE THE EQUATION OF THE DATA POINT
  195 W(6)=YD(IM)
      C=WD(IM)**2
      IM=IM+1
C     REVISE NC PRIOR TO THE NEXT GIVENS TRANSFORMATION
  200 NC=NC-1
      IF (NC.LE.0) GO TO 160
C     MAKE THE GIVENS TRANSFORMATION
  210 ALPHA=C*W(2)
      IF (ALPHA.NE.ZERO) GO TO 230
      DO 220 I=1,NC
  220 W(I+1)=W(I+2)
      GO TO 250
  230 DSTAR=W(KR)+ALPHA*W(2)
      BETA=W(KR)/DSTAR
      C=BETA*C
      GAMMA=ALPHA/DSTAR
      W(KR)=DSTAR
      DSTAR=W(2)
      DO 240 I=1,NC
      W(I+1)=W(I+2)-DSTAR*W(KR+I)
  240 W(KR+I)=BETA*W(KR+I)+GAMMA*W(I+2)
  250 KR=KR+7
      GO TO 200
C     ADVANCE THE CURRENT ROW OF THE UPPER TRIANGULAR MATRIX
  260 IF (IN.GE.NN) GO TO 330
      W(KK+28)=ZERO
      DO 270 I=4,28,6
      W(KK+I+1)=W(KK+I)
  270 W(KK+I)=ZERO
      KK=KK+7
C     CALCULATE THE THIRD DERIVATIVES OF THE B-SPLINES AT XN(IN+)
  280 W(2)=SIX/(W(IN+19)-W(IN+18))
      DO 300 J=2,4
      I=J+1
      W(I)=ZERO
  290 II=IN+I+17
      W(I)=(W(I-1)-W(I))/(W(II)-W(II-J))
      I=I-1
      IF (I.GE.2) GO TO 290
  300 CONTINUE
      IF (SWITCH) GO TO 130
C     SET UP THE EQUATION FOR THE SMOOTHING TERM AT A KNOT
      ALPHA=ZERO
      I=6
  310 W(I)=W(I-1)-ALPHA
      I=I-1
      IF (I.LE.1) GO TO 320
      ALPHA=W(I+6)
      W(I+6)=W(I)
      GO TO 310
  320 C=THETA(IN)**2
      IN=IN+1
      GO TO 210
C     BACK-SUBSTITUTE TO OBTAIN THE MULTIPLIERS OF THE B-SPLINES
  330 KU=KK
      KUU=KU+21
      KB=KS+7
      IBS=5
      DO 335 I=3,6
      KK=KK+7
  335 W(KK-2)=W(KK-I)
  340 KK=KUU+IBS
      KR=KK
  350 KK=KK-7
      K=KK
      J=KK-IBS
  360 J=J+1
      K=K+7
      W(KK)=W(KK)-W(J)*W(K)
      IF (K.LT.KR) GO TO 360
      KR=MIN0(KR,KK+21)
      IF (KK.GT.KB) GO TO 350
      IF (IBS-6) 370,470,500
  370 SWITCH=.TRUE.
      ALPHA=ONE
      DO 380 K=KS,KUU,7
  380 ALPHA=DMIN1(ALPHA,W(K))
      IF (ALPHA.LE.ZERO) GO TO 592
      IF (NN.LE.2) GO TO 600
C     NEXT A SPLINE WITH FIXED VALUES OF S'''(XN(1)) AND S'''(XN(N))
C     IS CALCULATED TO OBTAIN THE OPTIMAL EXTREME VALUES OF S'''(X).
C     SET THE COEFFICIENTS OF S'''(X(1)) AND S'''(X(N))
      KR=KS
      KK=KU
      DO 410 I=8,11
      W(KK+4)=W(I)
      KK=KK+7
      W(KR+6)=W(I+4)
  410 KR=KR+7
      DO 420 K=KR,KUU,7
  420 W(K+6)=ZERO
      IBS=6
      KRR=KS
C     BACK-SUBSTITUTE USING THE TRANSPOSE OF THE UPPER TRIANGULAR
C     MATRIX TO CALCULATE THE RESIDUALS OF EACH EXTRA SPLINE
  425 I=KRR+7
      KR=KRR
      DO 450 KK=I,KUU,7
      J=KK
      K=KK
  430 K=K-7
      J=J-6
      W(KK+IBS)=W(KK+IBS)-W(J)*W(K+IBS)
      IF (K.GT.KR) GO TO 430
  450 KR=MAX0(KR,KK-21)
      DO 460 KK=KRR,KUU,7
  460 W(KK+IBS)=W(KK+IBS)/W(KK)
      IF (IBS.EQ.6) GO TO 340
      DO 465 K=KS,KUU,7
      W(K+7)=ZERO
  465 IF (K.GE.KU) W(K+7)=W(K+4)
      IBS=7
      GO TO 340
C     BACK-SUBSTITUTE TO GIVE THE VECTOR FOR S'''(X(N))
  470 IBS=4
      KRR=KU
      GO TO 425
C     AT EACH DATA POINT CALCULATE THE RESIDUALS OF THREE SPLINES
  500 DO 520 I=6,10
  520 W(I)=ZERO
      IN=2
      KK=KS
      DO 580 IM=1,MM
  530 IF (XD(IM).LE.XN(IN)) THEN
        W(2)=ONE/(W(IN+18)-W(IN+17))
        DO 534 J=2,4
        I=J+1
        W(I)=ZERO
  532   II=IN+I+16
        C=(XD(IM)-W(II-J))*W(I-1)+(W(II)-XD(IM))*W(I)
        IF (DABS(C).LE.PREC) C=0.
        W(I)=C/(W(II)-W(II-J))
        I=I-1
        IF (I.GE.2) GO TO 532
  534   CONTINUE
        IF (SWITCH) GO TO 540
        GO TO 195
      ELSE
        IN=IN+1
        KK=KK+7
        GO TO 530
      ENDIF
  540 KR=KK
      W(11)=YD(IM)
      W(12)=ZERO
      W(13)=ZERO
      DO 560 I=2,5
      DO 550 J=5,7
  550 W(J+6)=W(J+6)-W(I)*W(KR+J)
  560 KR=KR+7
      DO 570 I=11,13
  570 W(I)=WD(IM)*W(I)
C     FORM THE NORMAL EQUATIONS TO FIX S'''(X) AT THE ENDS OF THE RANGE
      DO 580 I=12,13
      K=I+I-29
      DO 580 J=11,I
  580 W(K+J)=W(K+J)+W(I)*W(J)
C     CALCULATE THE B-SPLINE MULTIPLIERS OF THE REQUIRED FIT
      ALPHA=(W(8)*W(9)-W(6)*W(10))/(W(7)*W(10)-W(9)**2)
      BETA=(-W(8)-ALPHA*W(9))/W(10)
      DO 590 K=KS,KUU,7
  590 W(K+5)=W(K+5)+ALPHA*W(K+6)+BETA*W(K+7)
      GO TO 600
C     PRINT A DIAGNOSTIC MESSAGE IF THERE IS INSUFFICIENT DATA
  592 PRINT 594
  594 FORMAT (//5X,'THE RESULTS FROM VB06AD ARE UNRELIABLE BECAUSE ',
     1'THERE IS INSUFFICIENT DATA TO DEFINE THE SPLINE UNIQUELY')
C     CALCULATE THE REQUIRED VALUES OF S(X) AND S'(X) AT THE KNOTS
  600 KK=KS
      DO 620 IN=1,NN
      FN(IN)=ZERO
      GN(IN)=ZERO
      W(2)=(W(IN+19)-W(IN+18))/((W(IN+19)-W(IN+17))*(W(IN+19)-W(IN+16)))
      W(3)=(W(IN+18)-W(IN+17))/((W(IN+19)-W(IN+17))*(W(IN+20)-W(IN+17)))
      W(4)=ZERO
      DO 610 I=1,3
      J=IN+I+14
      ALPHA=((XN(IN)-W(J))*W(I)+(W(J+4)-XN(IN))*W(I+1))/(W(J+4)-W(J))
      FN(IN)=FN(IN)+ALPHA*W(KK+5)
      BETA=3.D0*(W(I)-W(I+1))/(W(J+4)-W(J))
      GN(IN)=GN(IN)+BETA*W(KK+5)
  610 KK=KK+7
  620 KK=KK-14
C     CALCULATE THE THIRD DERIVATIVE DISCONTINUITIES OF THE SPLINE
      DN(1)=ZERO
      IN=1
      IM=1
  630 H=XN(IN+1)-XN(IN)
      ALPHA=FN(IN)-FN(IN+1)+H*GN(IN)
      BETA=ALPHA+FN(IN)-FN(IN+1)+H*GN(IN+1)
      DN(IN+1)=SIX*BETA/H**3
      DN(IN)=DN(IN+1)-DN(IN)
      IN=IN+1
C     CALCULATE THE RESIDUALS AT THE DATA POINTS
  640 IF (IM.GT.MM) GO TO 650
      IF (XD(IM).GT.XN(IN)) GO TO 630
      C=(XD(IM)-XN(IN-1))/H
      RD(IM)=YD(IM)-C*FN(IN)+(C-ONE)*(FN(IN-1)+C*(ALPHA-C*BETA))
      IM=IM+1
      GO TO 640
  650 IF (IN.LT.NN) GO TO 630
C     PROVIDE PRINTING IF REQUESTED
      IF (IPRINT.LE.0) GO TO 770
      PRINT 660
  660 FORMAT (1H1,35X,'SPLINE APPROXIMATION OBTAINED BY VB06AD'//4X,'I',
     19X,'XN(I)',18X,'FN(I)',18X,'GN(I)',13X,'3RD DERIV CHANGE',
     211X,'THETA(I)'//)
      I=1
  670 PRINT 680,I,XN(I),FN(I),GN(I)
  680 FORMAT (I5,5D23.14)
  690 I=I+1
      IF (I-NN) 700,670,710
  700 PRINT 680,I,XN(I),FN(I),GN(I),DN(I),THETA(I)
      GO TO 690
  710 PRINT 720
  720 FORMAT (///4X,'I',9X,'XD(I)',18X,'YD(I)',18X,'WD(I)',19X,'FIT',
     118X,'RESIDUAL'//)
      IN=2
      DO 760 IM=1,MM
  730 IF (XD(IM).LE.XN(IN)) GO TO 750
      IN=IN+1
      IF (SWITCH) GO TO 730
      SWITCH=.TRUE.
      PRINT 740
  740 FORMAT (5X)
      GO TO 730
  750 C=YD(IM)-RD(IM)
      PRINT 680,IM,XD(IM),YD(IM),WD(IM),C,RD(IM)
  760 SWITCH=.FALSE.
  770 RETURN
      END
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS VC03AD
C###### CALLS   VB06
C  *********************************************************************
C  *                                                                   *
C  *  DUMMY: This file contains dummies for all the external routines  *
C  *  used throughout LIM3.  The program was designed and optimised    *
C  *  on the IBM3090E at JET and makes use of some IBM utilities,      *
C  *  the Harwell subroutine library ("HSL") and the NOCORONA package  *
C  *  also developed at JET.                                           *
C  *    The code collected here will require modifications when        *
C  *  transporting between machines, whilst the remainder of the LIM3  *
C  *  code should remain intact.                                       *
C  *                                                                   *
C  *  Routines in this file :-                                         *
C  *    XUFLOW: Prevent underflow interrupts                           *
C  *    SURAND: Get vector of random numbers in (0,1)                  *
C  *    ZA08AS: Get time of day as 8 character string                  *
C  *    ZA09AS: Get date as 8 character string                         *
C  *    ZV01AD: Get name of datafile connected to channel 5            *
C  *    TB04A : Create cubic spline                                    *
C  *    TG01B : Extract value along newly created cubic spline         *
C  *    VC03A : Produce best fit curve through data                    *
C  *    ZA02AS: Get time used so far in seconds                        *
C  *    RANINI: Dummy interface to RANSET random number initialiser    *
C  *                                                                   *
C  *                                      C.M.Farrell   June 1988      *
C  *                                                                   *
C  *********************************************************************
C
C     XUFLOW
C     ======
C     IBM  : System routine to prevent underflow interrupts occuring
C     CRAY : Replace with dummy routine here.
C
      SUBROUTINE XUFLOW (IFLAG)
      INTEGER IFLAG
      WRITE (6,'('' XUFLOW: dummied out for this application.'')')
      RETURN
      END
C
C***********************************************************************
C
C     ZA08AS
C     ======
C     IBM  : Harwell library routine to extract Time in 8 characters
C     CRAY : Replace with call to CLOCK system routine.
C
      SUBROUTINE ZA08AS (SYSTIM)
      CHARACTER*8 SYSTIM
       CALL CLOCKTEST (SYSTIM)
      RETURN
      END
C
C***********************************************************************
C
C     ZA09AS
C     ======
C     IBM  : Harwell library routine to extract Date in 8 characters
C     CRAY : Replace with call to DATE system routine.
C
      SUBROUTINE ZA09AS (SYSDAT)
      CHARACTER*8 SYSDAT
      CALL DATETEST (SYSDAT)
C      SYSDAT = '04/11/88'
      RETURN
      END
C
C***********************************************************************
C
C     ZV01AD
C     ======
C     IBM  : Harwell library routine to extract Dataset name connected
C            to channel 5, and to put the name in JFCB(1:44) and the
C            member in JFCB(45:52).
C     CRAY : No equivalent system routine.
C
      SUBROUTINE ZV01AD (IUNIT, VSN, DSN, JFCB)
      CHARACTER VSN*8,DSN(3)*8,JFCB*176
      WRITE (6,'('' ZV01AD: dummied out for this application.'')')
      JFCB = 'SYSTEM.TEST.FOR.IBM(TESTONE)      '
      RETURN
      END
C
C***********************************************************************
C
C     TB04A
C     =====
C     IBM  : Harwell library routine to create a cubic spline
C     CRAY : Use HSL routine with double precision name TB04AD.
C
      SUBROUTINE TB04A (N,X,F,FDASH,WORK)
      INTEGER N,I
      REAL X(N),F(N),FDASH(N),WORK(3*N)
C
      AUTOMATIC X2,F2,FDASH2,WORK2
      DOUBLE PRECISION X2(1001),F2(1001),FDASH2(1001),WORK2(3003)
C
C     IN ORDER TO INTERFACE THESE ROUTINES PROPERLY, SINCE
C     FORTRAN DOES NOT DO ARGUMENT TYPE CHECKING IT IS 
C     NECESSARY TO CONVERT THE SINGLE PRECISION VARIABLES TO 
C     DOUBLE PRECISION. THIS IS COMMENTED OUT FOR CRAY USAGE 
C     SINCE ONLY REALS ARE USED IN THIS ENVIRONMENT (64 BIT)
C
      DO 10 I = 1,N
         X2(I) = DBLE(X(I))
         F2(I) = DBLE(F(I))
         FDASH2(I) = DBLE(FDASH(I))
         WORK2(I) = DBLE(WORK(I))
         WORK2(N+I) = DBLE(WORK(N+I))
         WORK2(2*N+I) = DBLE(WORK(2*N+I))
10    CONTINUE         

      CALL TB04AD(N,X2,F2,FDASH2,WORK2)

C      CALL TB04AD (N,X,F,FDASH,WORK)

      DO 20 I = 1,N
         X(I) = SNGL(X2(I))
         F(I) = SNGL(F2(I))
         FDASH(I) = SNGL(FDASH2(I))
         WORK(I) = SNGL(WORK2(I))
         WORK(N+I) = SNGL(WORK2(N+I))
         WORK(2*N+I) = SNGL(WORK2(2*N+I))
20    CONTINUE         

      RETURN
      END
C
C***********************************************************************
C
C     TG01B
C     =====
C     IBM  : Harwell library function to extract value along spline.
C     CRAY : Use HSL function with double precision name TG01BD.
C
      REAL FUNCTION TG01B (IFLAG,N,X,F,FDASH,X0)
      REAL X(N),F(N),FDASH(N),X0
C     AUTOMATIC X2,F2,FDASH2,X02
      DOUBLE PRECISION X2(1001),F2(1001),FDASH2(1001),X02,TG01BD
      EXTERNAL TG01BD 
      INTEGER I
C
C     THIS IS AN INTERFACE BETWEEN A REAL AND DOUBLE PRECISION 
C     FUNCTION. THE CONVERSIONS WERE NOT REQUIRED ON A CRAY WHERE
C     A REAL WAS EQUIVALENT TO DOUBLE PRECISION ON A SMALLER MACHINE
C     AND THE CRAY DOUBLE PRECISION WAS DEFAULTED TO REAL. HOWEVER,
C     THIS DIFFERENCE IS SIGNIFICANT ON A WORKSTATION. 
C
C     D.ELDER   1991 APRIL 24
C
      X02 = dble(X0)
      DO 10 I = 1,N
         X2(I) = dble(X(I))
         F2(I) = dble(F(I))
         FDASH2(I) = dble(FDASH(I))
10    CONTINUE      

      TG01B = SNGL(TG01BD (IFLAG,N,X2,F2,FDASH2,X02))

      X0 = sngl(X02)
      DO 20 I = 1,N
         X(I) = sngl(X2(I))
         F(I) = sngl(F2(I))
         FDASH(I) = sngl(FDASH2(I))
20    CONTINUE      

      RETURN
      END
C
C***********************************************************************
C
C     VC03A
C     =====
C     IBM  : Harwell library routine to produce the curve of best-fit
C            through a collection of data points  (by a cubic spline
C            method).
C     CRAY : Use HSL routine with double precision name VC03AD.
C
      SUBROUTINE VC03A (M,N,XD,YD,WD,RD,XN,FN,GN,DN,THETA,IPRINT,W)
      INTEGER M,N,IPRINT,I,J,K
      REAL XD(M),YD(M),WD(M),RD(M),XN(N),FN(N),GN(N),DN(N),THETA(N)
      REAL W(7*M+9*N+6)
      AUTOMATIC  XD2,YD2,WD2,RD2,XN2,FN2
      AUTOMATIC GN2,DN2,THETA2,W2 
      DOUBLE PRECISION XD2(101),YD2(101),WD2(101),RD2(101),XN2(101)
      DOUBLE PRECISION FN2(101)
      DOUBLE PRECISION GN2(101),DN2(101),THETA2(101),W2(2002) 
C
C     THIS IS AN INTERFACE BETWEEN A REAL AND DOUBLE PRECISION 
C     FUNCTION. THE CONVERSIONS WERE NOT REQUIRED ON A CRAY WHERE
C     A REAL WAS EQUIVALENT TO DOUBLE PRECISION ON A SMALLER MACHINE
C     AND THE CRAY DOUBLE PRECISION WAS DEFAULTED TO REAL. HOWEVER,
C     THIS DIFFERENCE IS SIGNIFICANT ON A WORKSTATION. 
C
C     D.ELDER   1991 APRIL 24
C
      DO 10 I = 1,M
         XD2(I) = XD(I)
         YD2(I) = YD(I)
         WD2(I) = WD(I)
         RD2(I) = RD(I)
10    CONTINUE
      DO 20 I = 1,N
         XN2(I) = XN(I)
         FN2(I) = FN(I)
         GN2(I) = GN(I)
         DN2(I) = DN(I)
         THETA2(I)= THETA(I)
20    CONTINUE             
      DO 30 I = 1,(7*M+9*N+6)
         W2(I) = W(I)
30    CONTINUE

      CALL VC03AD (M,N,XD2,YD2,WD2,RD2,XN2,FN2,GN2,DN2,
     >             THETA2,IPRINT,W2)

      DO 40 I = 1,M
         XD(I) = XD2(I)
         YD(I) = YD2(I)
         WD(I) = WD2(I)
         RD(I) = RD2(I)
40    CONTINUE
      DO 50 I = 1,N
         XN(I) = XN2(I)
         FN(I) = FN2(I)
         GN(I) = GN2(I)
         DN(I) = DN2(I)
         THETA2(I)= THETA(I)
50    CONTINUE             
      DO 60 I = 1,(7*M+9*N+6)
         W(I) = W2(I)
60    CONTINUE

      RETURN
      END
C
C***********************************************************************
C
C     ZA02AS
C     ======
C     IBM  : Harwell library function to extract CPU time used so far.
C     HOT  : For hotspot analysis, dummy out by setting ZA02AS = 0.0
C     CRAY : Replace with system function SECOND.
C
      REAL FUNCTION ZA02AS (IFLAG)
      INTEGER I
      I = MCLOCK()
      ZA02AS = I/100.0
CHOT  ZA02AS = 0.0
C     ZA02AS = SECOND ()
      RETURN
      END
C
C***********************************************************************
C
C     RANINI
C     ======
C     IBM  : Dummy routine - no need to call RANSET
C     CRAY : Interface to random no. initialiser system routine RANSET
C
      SUBROUTINE RANINI (ISEED)
      INTEGER ISEED
      CALL NEWSRAND (ISEED)
C
C     CALL SRAND (ISEED ) FOR CRAY OR IBM DEFAULT GENERATOR
C
      RETURN
      END
C
C***********************************************************************
C
C     SURAND
C     ======
C     IBM  : ESSL library routine to generate vector of random numbers
C     CRAY : Replace with calls to RANF generator within a vectorisable
C            loop.
C
      SUBROUTINE SURAND (SEED,NRANDS,RANDS)      
      DOUBLE PRECISION SEED
      INTEGER NRANDS  
      REAL RANDS(NRANDS),NEWRAND
      EXTERNAL NEWRAND
      DO 100 J = 1, NRANDS
          RANDS(J) = NEWRAND()
C         
C         RANJ(J) = RAND ()  FOR CRAY OR IBM USING DEFAULT GENERATOR
C
  100 CONTINUE
      RETURN
      END

