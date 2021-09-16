C     -*-Fortran-*-
C
C***********************************************************************
C
C     TB04A
C     =====
C     IBM  : HARWELL LIBRARY ROUTINE TO CREATE A CUBIC SPLINE
C     CRAY : USE HSL ROUTINE WITH DOUBLE PRECISION NAME TB04AD.
C
      SUBROUTINE TB04A (N,X,F,FDASH,WORK)
      implicit none
      INTEGER N,I
      REAL X(N),F(N),FDASH(N),WORK(3*N)
C
C     AUTOMATIC X2,F2,FDASH2,WORK2
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
C     IBM  : HARWELL LIBRARY FUNCTION TO EXTRACT VALUE ALONG SPLINE.
C     CRAY : USE HSL FUNCTION WITH DOUBLE PRECISION NAME TG01BD.
C
      REAL FUNCTION TG01B (IFLAG,N,X,F,FDASH,X0)
      implicit none
      integer n,iflag
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
      X02 = DBLE(X0)
      DO 10 I = 1,N
         X2(I) = DBLE(X(I))
         F2(I) = DBLE(F(I))
         FDASH2(I) = DBLE(FDASH(I))
10    CONTINUE

      TG01B = SNGL(TG01BD (IFLAG,N,X2,F2,FDASH2,X02))

      X0 = SNGL(X02)
      DO 20 I = 1,N
         X(I) = SNGL(X2(I))
         F(I) = SNGL(F2(I))
         FDASH(I) = SNGL(FDASH2(I))
20    CONTINUE

      RETURN
      END
C
C***********************************************************************
C
C     VC03A
C     =====
C     IBM  : HARWELL LIBRARY ROUTINE TO PRODUCE THE CURVE OF BEST-FIT
C            THROUGH A COLLECTION OF DATA POINTS  (BY A CUBIC SPLINE
C            METHOD).
C     CRAY : USE HSL ROUTINE WITH DOUBLE PRECISION NAME VC03AD.
C
      SUBROUTINE VC03A (M,N,XD,YD,WD,RD,XN,FN,GN,DN,THETA,IPRINT,W)
      implicit none
      INTEGER M,N,IPRINT,I,J,K
      REAL XD(M),YD(M),WD(M),RD(M),XN(N),FN(N),GN(N),DN(N),THETA(N)
      REAL W(7*M+9*N+6)
C     AUTOMATIC  XD2,YD2,WD2,RD2,XN2,FN2
C     AUTOMATIC GN2,DN2,THETA2,W2
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
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS TB04AD
      SUBROUTINE TB04AD(N,X,F,D,A)
      implicit none
      integer n
C STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
      DOUBLE PRECISION A,D,F,H1,H2,P,X
      DIMENSION X(N),F(N),D(N),A(*)
      integer np
      DATA NP/6/
c
c     Local variables
c
      integer i,j,k
c
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
      implicit none
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
      integer iflg,ieps,ix,n,j
      DOUBLE PRECISION A,B,D,H,Q1,Q2,S,SS,U,X,Z
C
C     ALLOWABLE ROUNDING ERROR ON POINTS AT EXTREAMS OF KNOT RANGE
C     IS 2**IEPS*MAX(!U(1)!,!U(N)!).
c slmod begin
c...  How did this ever work? -SL 19.08.06
      DIMENSION U(N),S(N),D(N)
c
c      DIMENSION U(1),S(1),D(1)
c slmod end
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
      implicit none  
      integer m,n,iprint
      DOUBLE PRECISION DN,FN,GN,HLF,HS,PRP,RD,RP,SA,SR,SW,THETA,
     *                 TMAX,W,WD,XD,XN,YD
      DIMENSION XD(1),YD(1),WD(1),RD(1),XN(*),FN(*),GN(1),DN(1),
     1THETA(1),W(1)
c
c     Local variables
c
      integer ndimax,mm,i,j,k,ja,iw,ip,nn,jj,iis,kc,jw,ik
      integer kl,ku,kku,kkl,kz
c
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
      W(J+2)=DBLE(I)
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
      RP=DSQRT(SR/DBLE(KC))
      IF (DBLE(KC)-SW) 32,32,34
   34 GO TO (35,36,37),IIS
   35 PRP=RP
      IIS=3
      IF (DBLE(KC)-2.D0*SW) 38,38,39
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
      implicit none 
      integer m,n,iprint
      DOUBLE PRECISION ALPHA,BETA,C,DN,DSTAR,FN,GAMMA,GN,H,PREC
      DOUBLE PRECISION RD,THETA,W,WD,XD,XN,YD,ZERO,ONE,TWO,SIX
      LOGICAL SWITCH
      DIMENSION XD(1),YD(1),WD(1),RD(1),XN(1),FN(1),GN(1),DN(1),
     1THETA(1),W(*)
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,SIX/6.0D0/
c
c     Local variables
c
      integer i,nn,mm,in,ks,ku,im,kk,kr,nc,j,ii,kuu,kb,ibs,k,krr
c
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
C######DATE   12 DEC 1990,     COPYRIGHT UKAEA, HARWELL.
C######ALIAS GA15A
C######CALLS FD05 KB05 NB01
C######DATE   24 JUN 1985     COPYRIGHT UKAEA, HARWELL.
C######ALIAS GA15A
C###### CALLS KB05  NB01
      SUBROUTINE GA15A(N, KIND, W, IW, IND, LND, X, Y, TT, XD, YD, LP)
      implicit none  
      integer n,iw,ind,lnd,lp,kind
      real w,x,y,xd,yd,tt  
      DIMENSION W(IW), X(N), Y(N), XD(N), YD(N), TT(N)
      DIMENSION IND(2,LND)
c
c     Local variables
c 
      integer i,i1,m,i2,i3,i4,n1,nstat,il
      real cubic,r,d,c,b,a,eps,fd05a,big,dy,dyn,t1,t2,xmax,xmin,f,dm
      integer ir,nrep,nnd,ncrit,ist,inc,k,j,in,in1
      real dp,root,root1,rt,yl,tl,yn,tn,yp,t,err,xn,xl
c
      real zero,one,two,three,four,eight
      DATA ZERO /0./, ONE /1./,  TWO /2.0/, THREE /3.0/,
     * FOUR /4.0/,  EIGHT /8.0/
c
      CUBIC(R) = D + R*(C+R*(B+R*A))
C  CHECK N GREATER THAN 4
C     write(6,*) 'GA15a:Called:', n,kind,iw,lnd,lp
      EPS=FD05A(1)*10.
      BIG=FD05A(5)*.1
      IF (N.LT.4) GO TO 430
C  CHECK DATA IS PERIODIC
      IF (X(1).NE.X(N) .OR. Y(1).NE.Y(N)) GO TO 420
C  CHECK KIND=1 OR KIND=2
      IF (KIND.NE.1 .AND. KIND.NE.2) GO TO 440
      IF (KIND.EQ.1) GO TO 20
C  CHECK VALUES OF TT IN ORDER (CUBICS ONLY)
      DO 10 I=2,N
        IF (TT(I).LE.TT(I-1)) GO TO 450
   10 CONTINUE
C  CHECK XD AND YD ARE PERIODIC  (CUBICS ONLY)
      IF (XD(1).NE.XD(N) .OR. YD(1).NE.YD(N)) GO TO 460
C
C     PARTITION THE STORAGE IN W
C     IN THE PARAMETRIC CUBIC CASE W(I),W(N+I),I=2,N HOLD THE MINIMUM
C  AND MAXIMUM VALUES OF X IN THE INTERVAL BETWEEN KNOTS I-1 AND I.
   20 I1 = 0
      IF (KIND.NE.1) I1 = 2*N
C     W(I1+I),I=1,2,... HOLD THE STATIONARY VALUES OF Y IN ASCENDING
C   ORDER.
      M = (IW-I1)/4
      I2 = I1 + M
C     W(I2+I),I=1,2,... HOLD THE VALUES OF CAPITAL T AT CRITICAL
C  POINTS IN SEQUENCE ROUND THE CURVE.
      I3 = I2 + M
      I4 = I3 + M
C     W(I3+I),W(I4+I),I=1,2,... HOLD THE MAXIMUM AND MINIMUM X VALUES
C  IN THE BOUNDARY SEGMENT ENDING AT THE ITH CRITICAL POINT.
      N1 = N - 1
      NSTAT = 0
      IF (KIND.NE.1) GO TO 50
C
C     FIND MAXIMA AND MINIMA OF Y IN POLYGONAL CASE
C     THE MAXIMA AND MINIMA OF Y AND CORRESPONDING CAPITAL T VALUES
C  ARE STORED TEMPORARILY IN W(I3+I),W(I4+I),I=1,2,...
      DY = Y(N) - Y(N-1)
      DO 40 I=2,N
        DYN = Y(I) - Y(I-1)
        IF (DY*DYN.GT.ZERO) GO TO 30
        NSTAT = NSTAT + 1
        IF (I3+NSTAT.GT.I4) GO TO 400
        W(I3+NSTAT) = Y(I-1)
        W(I4+NSTAT) = I - 1
   30   DY = DYN
   40 CONTINUE
      GO TO 140
C
C     FIND MAXIMA AND MINIMA OF X AND Y.
C     THE MAXIMA AND MINIMA OF Y AND CORRESPONDING CAPITAL T VALUES
C  ARE STORED TEMPORARILY IN W(I3+I),W(I4+I),I=1,2,...
   50 DO 130 I=2,N
        T1 = I - 1
        T2 = I
        CALL GA15C(T1, T2, X, TT, XD, XMAX, XMIN, N)
        W(N+I) = XMAX
        W(I) = XMIN
        F = (Y(I)-Y(I-1))/FOUR
        A = (TT(I)-TT(I-1))/EIGHT
        DM = YD(I-1)*A
        DP = YD(I)*A
        A = THREE*(DM+DP-F)
        B = DP - DM
        C = THREE*F - DM - DP
        ROOT = -TWO
        ROOT1 = TWO
        IF (B*B.LE.A*C) GO TO 60
        D = B + SIGN(SQRT(B*B-A*C),B)
        ROOT = -C/D
        IF (A.NE.ZERO) ROOT1 = -D/A
C     TEST FOR STATIONARY POINTS AT ENDS OF RANGE.
   60   D = (ABS(A)+TWO*ABS(B)+ABS(C))*EPS
        IF (FOUR*ABS(DP).GT.D) GO TO 80
        IF (ABS(ROOT-ONE).LT.ABS(ROOT1-ONE)) GO TO 70
        ROOT1 = ONE
        GO TO 80
   70   ROOT = ONE
   80   IF (FOUR*ABS(DM).GT.D) GO TO 100
        IF (ABS(ROOT+ONE).LT.ABS(ROOT1+ONE)) GO TO 90
        ROOT1 = -ONE
        GO TO 100
   90   ROOT = -ONE
C     STORE STATIONARY POINTS
  100   A = DM + DP - F
        D = Y(I) - A - B - C
        RT =AMIN1(ROOT,ROOT1)
        DO 120 IR=1,2
          IF (ABS(RT).GT.ONE) GO TO 110
          NSTAT = NSTAT + 1
          IF (NSTAT+I3.GT.I4) GO TO 400
          W(NSTAT+I4) = I + (RT-ONE)/TWO
          W(NSTAT+I3) = CUBIC(RT)
  110     RT = AMAX1(ROOT,ROOT1)
  120   CONTINUE
  130 CONTINUE
C
C     SORT STATIONARY POINTS INTO ASCENDING ORDER.
  140 IND(2,1) = NSTAT
      IF (NSTAT.GT.LND) GO TO 410
      DO 150 I=1,NSTAT
        IND(2,I) = -1
        W(I+I1) = W(I+I3)
  150 CONTINUE
      CALL KB05A(W(I1+1), NSTAT)
C
C     REMOVE ANY REPEATED STATIONARY POINTS.
      NREP = 0
      DO 160 I=2,NSTAT
        IF (W(I1+I-1).EQ.W(I1+I)) NREP = NREP + 1
        W(I1+I-NREP) = W(I1+I)
  160 CONTINUE
C
C     FIND CRITICAL POINTS AND INDEX THE BOUNDARY INTERVALS.
C     FOR EACH INTERVAL BETWEEN W(I1+I-1) AND W(I1+I),I=2,3,...
C  AN ASSOCIATED LIST OF INTEGERS IS HELD IN IND STARTING AT IND(1,I).
C  THE INTEGERS THEMSELVES ARE IN IND(1,.) AND WITH EACH IS
C  ASSOCIATED THE ADDRESS, HELD IN IND(2,.) OF THE NEXT INTEGER OR
C  ZERO TO INDICATE THE END.  EACH INTEGER K IN THIS LIST POINTS
C  TO AN INTERVAL W(I2+IABS(K)-1) TO W(I2+IABS(K)), WHOSE Y VALUES SPAN
C  THE INTERVAL W(I1+I-1) TO W(I1+I).  Y IS INCREASING OR DECREASING
C  IN THIS INTERVAL ACCORDING AS K IS POSITIVE OR NEGATIVE.
      YL = W(I3+NSTAT)
      TL = W(I4+NSTAT)
      NND = NSTAT - NREP
      NCRIT = 0
      IST = 0
  170 IST = IST + 1
      IF (W(I1+IST).NE.YL) GO TO 170
      DO 310 I=1,NSTAT
        YN = W(I3+I)
        TN = W(I4+I)
        IF (TL.GT.TN) TN = TN + N1
        INC = 1
        IF (YN.LT.YL) INC = -1
  180   NCRIT = NCRIT + 1
        IF (I2+NCRIT.GT.I3) GO TO 400
        W(I2+NCRIT) = W(I4+I)
        IF (YN.EQ.YL) GO TO 300
        IST = IST + INC
        YP = W(I1+IST)
        IF (YP.EQ.YN) GO TO 250
        K = TL
  190   J = K
        K = K + 1
        IF (J.GE.N) J = J - N1
        IF ((Y(J+1)-YP)*INC.LE.ZERO .AND. K.LT.TN) GO TO 190
        T1 = AMAX1(J+ZERO,TL)
        T2 = AMIN1(J+ONE,TN)
        IF (KIND.EQ.1) GO TO 240
C     FIND THE INTERCEPT OF Y=YP WITH THE  SPLINE IN Y
C    IN THE INTERVAL T1,T2 (ASSUMED TO LIE IN ONE KNOT INTERVAL).
        J = T1
        A = (TT(J+1)-TT(J))/EIGHT
        DM = YD(J)*A
        DP = YD(J+1)*A
        F = (Y(J+1)-Y(J))/FOUR
        A = DM + DP - F
        B = DP - DM
        C = THREE*F - DM - DP
        D = Y(J+1) - A - B - C
        K = 0
        T = AMAX1(ABS(T1),ABS(T2))*TWO
        ERR = EPS*(ABS(D)+ABS(C)*(ONE+T)+ABS(B)*(ONE+T+T)+ABS(A)*
     *   (ONE+T+T+T))
  200   CALL NB01A(K, T1, T2, ERR, T, R, 100)
        GO TO (210, 230, 220, 220), K
  210   R = CUBIC((T-J)*TWO-ONE) - YP
        GO TO 200
  220   WRITE (LP,99999)
        IND(1,1) = 4
        GO TO 390
  230   W(I2+NCRIT) = T
        GO TO 250
  240   W(I2+NCRIT) = J + 1 + (YP-Y(J+1))/(Y(J+1)-Y(J))
C     INSERT CRITICAL INTERVAL IN APPROPRIATE LIST.
  250   J = IST - (INC-1)/2
  260   IF (IND(2,J)) 290, 280, 270
  270   J = IND(2,J)
        GO TO 260
  280   IND(2,J) = NND + 1
        NND = NND + 1
        IF (NND.GT.LND) GO TO 410
        J = NND
  290   IND(1,J) = NCRIT*INC
        IND(2,J) = 0
        IF (YP.NE.YN) GO TO 180
  300   TL = W(I4+I)
        YL = YN
  310 CONTINUE
      NSTAT = NSTAT - NREP
C
C     FIND MAXIMUM AND MINIMUM X VALUES FOR EACH INTERVAL BETWEEN
C  CRITICAL POINTS.
      TL = W(I2+NCRIT)
      DO 380 I=1,NCRIT
        TN = W(I2+I)
        IN = TN
        IN1 = IN + 1
        IF (IN1.GT.N) IN1 = 2
        XN = X(IN) + (TN-IN)*(X(IN1)-X(IN))
        IF (TL.GT.TN) TN = TN + N1
        IL = TL
        IN = TN
        IF (KIND.NE.1) GO TO 330
        XL = X(IL) + (TL-IL)*(X(IL+1)-X(IL))
        W(I3+I) = AMAX1(XL,XN)
        W(I4+I) = AMIN1(XL,XN)
        IL = IL + 1
        IF (IL.GT.IN) GO TO 370
        DO 320 K=IL,IN
          J = K
          IF (J.GE.N) J = J - N1
          W(I3+I) = AMAX1(W(I3+I),X(J))
          W(I4+I) = AMIN1(W(I4+I),X(J))
  320   CONTINUE
        GO TO 370
  330   DO 360 J=IL,IN
          T1 = AMAX1(J+ZERO,TL)
          T2 = AMIN1(J+ONE,TN)
          IF (T2.LE.N) GO TO 340
          T1 = T1 - N1
          T2 = T2 - N1
  340     CALL GA15C(T1, T2, X, TT, XD, XMAX, XMIN, N)
          IF (J.NE.IL) GO TO 350
          W(I3+I) = XMAX
          W(I4+I) = XMIN
  350     W(I3+I) = AMAX1(W(I3+I),XMAX)
          W(I4+I) = AMIN1(W(I4+I),XMIN)
  360   CONTINUE
  370   TL = W(I2+I)
  380 CONTINUE
      IND(1,1) = -NSTAT
      IND(2,1) = NCRIT
      RETURN
C
C     ERROR RETURNS
  390 WRITE (LP,99998)
      RETURN
  400 WRITE (LP,99997)
      IND(1,1) = 1
      GO TO 390
  410 WRITE (LP,99996)
      IND(1,1) = 2
      GO TO 390
  420 WRITE (LP,99995)
      IND(1,1) = 3
      GO TO 390
  430 WRITE (LP,99994)
      IND(1,1) = 5
      GO TO 390
  440 WRITE (LP,99993)
      IND(1,1) = 6
      GO TO 390
  450 WRITE (LP,99992)
      IND(1,1) = 7
      GO TO 390
  460 WRITE (LP,99991)
      IND(1,1) = 8
      GO TO 390
99999 FORMAT (32X, 15HOF NB01 FAILURE)
99998 FORMAT (31H+ERROR RETURN FROM GA15 BECAUSE)
99997 FORMAT (32X, 15HIW IS TOO SMALL)
99996 FORMAT (32X, 16HLND IS TOO SMALL)
99995 FORMAT (32X, 24HTHE DATA IS NOT PERIODIC)
99994 FORMAT (32X, 14HN IS TOO SMALL)
99993 FORMAT (32X, 24HKIND NOT EQUAL TO 1 OR 2)
99992 FORMAT (32X, 14HT NOT IN ORDER)
99991 FORMAT (32X, 24HXD OR YD IS NOT PERIODIC)
      END
      SUBROUTINE GA15B(XP, YP, RESULT, N, KIND, W, IW, IND, LND, X, Y,
     * TT, XD, YD, LP)
      implicit none
c
      integer n,iw,ind,lnd,lp,kind
      real w,x,y,xd,yd,tt,result,xp,yp  
c
      DIMENSION W(IW), X(N), Y(N), XD(N), YD(N), TT(N)
      DIMENSION IND(2,LND)
c
c     Local variables
c 
c
c     Local variables
c 
      integer i,m,i2,i3,i4,n1,nstat
      real cubic,r,d,c,b,a,eps,fd05a,big,dy,dyn,t1,t2,xmax,xmin,f,dm
      integer ir,nrep,nnd,ncrit,ist,inc,k,j,in,in1,i1
      real dp,root,root1,rt,yl,tl,yn,tn,t,err,xn,xl
c
      integer jl,ju,jt,nl,iu,il,it
      real dist,xdist,resl
c
      real zero,one,two,three,four,eight
      DATA ZERO /0.0/, ONE /1.0/,  TWO /2.0/, THREE /3.0/,
     * FOUR /4.0/,  EIGHT /8.0/
      CUBIC(R) = D + R*(C+R*(B+R*A))
      BIG=FD05A(5)*.1
      EPS=FD05A(1)*10.
C     PARTITION THE STORAGE IN W
      I1 = 0
      IF (KIND.NE.1) I1 = 2*N
      M = (IW-I1)/4
      I2 = I1 + M
      I3 = I2 + M
      I4 = I3 + M
      N1 = N - 1
      NSTAT = -IND(1,1)
      NCRIT = IND(2,1)
C
C     FIND THE INTERVAL BETWEEN STATIONARY POINTS.
C     W(I1+JL).LT.YP.LE.W(I1+JU)
      JL = 0
      JU = NSTAT + 1
   10 IF (JU.EQ.JL+1) GO TO 30
      JT = (JL+JU)/2
      IF (YP.GT.W(I1+JT)) GO TO 20
      JU = JT
      GO TO 10
   20 JL = JT
      GO TO 10
   30 DIST = BIG
      RESULT = -BIG
      IF (JU.EQ.1) GO TO 220
      IF (JU.GT.NSTAT) GO TO 230
C
C     COUNT THE BOUNDARY SEGMENTS TO THE LEFT IN NL.  XDIST IS SET
C  POSITIVE OR NEGATIVE ACCORDING AS THE BOUNDARY SEGMENT IS TO THE LEFT
C  OR RIGHT.  ABS(XDIST) IS A LOWER BOUND ON DISTANCE TO THE BOUNDARY
C  AND IS ZERO FOR POINTS ON THE BOUNDARY.
   40 NL = 0
   50 J = IND(1,JU)
      F = J
      J = IABS(J)
C     TEST AGAINST MAXIMUM AND MINIMUM X VALUES FOR THE BOUNDARY
C  SEGMENT.
      XDIST = XP - W(I3+J)
      IF (XDIST.GT.ZERO) GO TO 200
      XDIST = XP - W(I4+J)
      IF (XDIST.LT.ZERO) GO TO 210
C     PERFORM A BINARY SEARCH TO FIND KNOT INTERVAL INVOLVED.
      IU = W(I2+J) + ONE
      IF (IU.EQ.N+1) IU = N
      IF (J.EQ.1) GO TO 60
      IL = W(I2+J-1)
      GO TO 70
   60 IL = W(I2+NCRIT)
   70 IF (IL.GE.IU) IL = IL - N1
      IF (IU.EQ.IL+1) GO TO 90
      IT = (IL+IU)/2
      IF (IT.LE.1) IT = IT + N1
      IF ((YP-Y(IT))*F.GT.ZERO) GO TO 80
      IU = IT
      IL = IL + N1
      GO TO 70
   80 IL = IT
      GO TO 70
C  CHECK FOR VERTEX
   90 XDIST=ZERO
      IF(XP.EQ.X(IU).AND.YP.EQ.Y(IU)) GO TO 190
      IF (KIND.NE.1) GO TO 100
C     PERFORM DIRECT TEST IN THE POLYGONAL CASE.
      XDIST = XP - (YP-Y(IL))*(X(IU)-X(IL))/(Y(IU)-Y(IL)) - X(IL)
      GO TO 190
C     TEST AGAINST MAXIMUM AND MINIMUM X VALUES IN KNOT INTERVAL.
  100 XDIST = XP - W(N+IU)
      IF (XDIST.GT.ZERO) GO TO 200
      XDIST = XP - W(IU)
      IF (XDIST.LT.ZERO) GO TO 210
C     FIND ACTUAL POINT OF INTERSECTION WITH THE BOUNDARY.
      T = W(I2+NCRIT)
      IF (J.NE.1) T = W(I2+J-1)
      T1 = IL
      IF (T.GT.IL .AND. T.LT.IU) GO TO 110
      XDIST = XP - X(IL)
      IF (YP-Y(IL)) 120, 190, 120
  110 T1 = T
  120 T = W(I2+J)
      T2 = IU
      IF (T.GT.IL .AND. T.LT.IU) GO TO 130
      XDIST = XP - X(IU)
      IF (YP.EQ.Y(IU)) GO TO 190
      GO TO 140
  130 T2 = T
C
C     FIND THE INTERCEPT OF Y=YP WITH THE  SPLINE IN Y
C    IN THE INTERVAL T1,T2 (ASSUMED TO LIE IN ONE KNOT INTERVAL).
  140 J = T1
      A = (TT(J+1)-TT(J))/EIGHT
      DM = YD(J)*A
      DP = YD(J+1)*A
      F = (Y(J+1)-Y(J))/FOUR
      A = DM + DP - F
      B = DP - DM
      C = THREE*F - DM - DP
      D = Y(J+1) - A - B - C
      K = 0
      T = AMAX1(ABS(T1),ABS(T2))*TWO
      ERR = EPS*(ABS(D)+ABS(C)*(ONE+T)+ABS(B)*(ONE+T+T)+ABS(A)*
     * (ONE+T+T+T))
  150 CALL NB01A(K, T1, T2, ERR, T, R, 100)
      GO TO (160, 180, 170, 170), K
  160 R = CUBIC((T-J)*TWO-ONE) - YP
      GO TO 150
  170 WRITE (LP,99999)
      IND(1,1) = 4
      WRITE (LP,99998)
      RETURN
  180 A = (TT(J+1)-TT(J))/EIGHT
      DM = XD(J)*A
      DP = XD(J+1)*A
      F = (X(J+1)-X(J))/FOUR
      A = DM + DP - F
      B = DP - DM
      C = THREE*F - DM - DP
      D = X(J+1) - A - B - C
      XDIST = XP - CUBIC((T-J)*TWO-ONE)
C     TEST XDIST
  190 IF (XDIST) 210, 240, 200
  200 NL = NL + 1
  210 DIST = AMIN1(DIST,ABS(XDIST))
      JU = IND(2,JU)
      IF (JU.NE.0) GO TO 50
C
C     RETURN OR HANDLE THE SPECIAL CASE WHERE YP EQUALS A STATIONARY Y
C  VALUE.
      RESULT = DIST
      IF ((NL/2)*2.EQ.NL) RESULT = -DIST
      IF (YP.EQ.W(I1+JL)) GO TO 250
      JU = JL + 1
  220 IF (YP.EQ.W(I1+JU)) GO TO 260
  230 RETURN
  240 RESULT = ZERO
      RETURN
  250 IF (RESULT*RESL) 240, 230, 230
  260 RESL = -ONE
      IF (JU.EQ.NSTAT) GO TO 250
      JL = JU
      JU = JL + 1
      IF (RESULT.GT.ZERO) RESL = ONE
      GO TO 40
99999 FORMAT (32X, 15HOF NB01 FAILURE)
99998 FORMAT (31H+ERROR RETURN FROM GA15 BECAUSE)
      END
      SUBROUTINE GA15C(T1, T2, X, TT, XD, XMAX, XMIN, N)
      implicit none 
      integer n
      real x,tt,xd,xmax,xmin
      DIMENSION X(N), XD(N), TT(N)
c
c     Local variables
c 
      integer l,ir
      real t1,t2,cubic,r,d,c,b,a,f,dm,dp,r1,r2,x1,x2,e,root,xx

c
      real zero,one,two,three,four,eight
      DATA ZERO /0.0/, ONE /1.0/, TWO /2.0/, THREE /3.0/, FOUR /4.0/,
     * EIGHT /8.0/
c
      CUBIC(R) = D + R*(C+R*(B+R*A))
C
C     THE FOLLOWING INSTRUCTIONS FIND THE MAXIMUM AND MINIMUM VALUES OF
C  X IN THE INTERVAL T1,T2 (ASSUMED TO LIE IN ONE KNOT INTERVAL).
      L = MIN0(N,INT(T1)+1)
      F = (X(L)-X(L-1))/FOUR
      A = (TT(L)-TT(L-1))/EIGHT
      DM = XD(L-1)*A
      DP = XD(L)*A
      A = DM + DP - F
      B = DP - DM
      C = THREE*F - DM - DP
      D = X(L) - A - B - C
      R1 = (T1-L)*TWO + ONE
      R2 = (T2-L)*TWO + ONE
      X1 = CUBIC(R1)
      X2 = CUBIC(R2)
      IF (R1.EQ.-ONE) X1 = X(L-1)
      IF (R2.EQ.ONE) X2 = X(L)
      XMAX = AMAX1(X1,X2)
      XMIN = AMIN1(X1,X2)
      IF (B*B.LE.THREE*A*C) GO TO 30
      E = B + SIGN(SQRT(B*B-THREE*A*C),B)
      ROOT = -C/E
      DO 20 IR=1,2
        IF (ROOT.LE.R1 .OR. ROOT.GE.R2) GO TO 10
        XX = CUBIC(ROOT)
        XMAX = AMAX1(XMAX,XX)
        XMIN = AMIN1(XMIN,XX)
   10   IF (A.NE.ZERO) ROOT = -E/A
   20 CONTINUE
   30 RETURN
      END
C######DATE   14 MAY 1991,     COPYRIGHT UKAEA, HARWELL.
C######ALIAS FD05A
*######DATE   27 FEB 1989     COPYRIGHT UKAEA, HARWELL.
*######ALIAS FD05A
      REAL FUNCTION FD05A( INUM )
      implicit none
      INTEGER INUM
      REAL RC( 5 )
C
C  REAL CONSTANTS (SINGLE PRECISION ARITHMETIC).
C
C  OBTAINED FROM H.S.L. SUBROUTINE ZE02AM.
C  NICK GOULD AND SID MARLOW, HARWELL, JULY 1988.
C
C  RC(1) THE SMALLEST POSITIVE NUMBER: 1.0 + RC(1) > 1.0.
C  RC(2) THE SMALLEST POSITIVE NUMBER: 1.0 - RC(2) < 1.0.
C  RC(3) THE SMALLEST NONZERO +VE REAL NUMBER.
C  RC(4) THE SMALLEST FULL PRECISION +VE REAL NUMBER.
C  RC(5) THE LARGEST FINITE +VE REAL NUMBER.
C
      DATA RC( 1 ) /    0.953674317E-06 /
      DATA RC( 2 ) /    0.596046449E-07 /
      DATA RC( 3 ) /    0.539760536E-37 /
      DATA RC( 4 ) /    0.539760536E-37 /
c
c      DATA RC( 3 ) /    0.539760536E-78 /
c      DATA RC( 4 ) /    0.539760536E-78 /
c
c     Patch for out of range error on SUN
c
      DATA RC( 5 ) /    0.723700514E+38 /
c
c      DATA RC( 5 ) /    0.723700514E+76 /
c
      IF ( INUM .LE. 0 .OR. INUM .GE. 6 ) THEN
         PRINT 2000, INUM
         STOP
      ELSE
         FD05A = RC( INUM )
      ENDIF
      RETURN
 2000 FORMAT( ' INUM =', I3, ' OUT OF RANGE IN FD05A.',
     *        ' EXECUTION TERMINATED.' )
      END
C######DATE   12 DEC 1990,     COPYRIGHT UKAEA, HARWELL.
C######ALIAS KB05A
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS KB05A
      SUBROUTINE KB05A(COUNT,N)
      implicit none
      integer n,mark
      real count
C
C             KB05A      HANDLES REAL SINGLE-LENGTH VARIABLES
C  THE WORK-SPACE 'MARK' OF LENGTH 50 PERMITS UP TO 2**(50/2) NUMBERS
C  TO BE SORTED. THIS IS MORE THAN THE IBM VIRTUAL MEMORY SPACE
C  WILL HOLD .
c slmod begin
c ARRAY BOUNDS: Generated an error since COUNT is indexed as an array below.
      DIMENSION COUNT(*),MARK(50)
c
c     Local variables:
c
      integer m,la,is,if,mloop,ifka,isl,j,iy,k,ifk,k1,ip
      integer lngth,ifend,ik,is1,i
      real av,x 
c
c      WRITE(0,*) 'MOD: MODIFICATION IN KB05A'
c
c      DIMENSION COUNT(1),MARK(50)
c slmod end
C  CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED
      IF(N.EQ.1)GOTO 280
      IF(N.GE.1)GO TO 110
      WRITE(6,100)
  100 FORMAT(///20X,65H ***KB05A*** NO NUMBERS TO BE SORTED ** RETURN TO
     2 CALLING PROGRAM )
      GOTO 280
C  'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER
C  THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.
  110 M=12
C  SET UP INITIAL VALUES.
      LA=2
      IS=1
      IF=N
      DO 270 MLOOP=1,N
C  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE .
      IFKA=IF-IS
      IF((IFKA+1).GT.M)GOTO 140
C********* FINAL SORTING ***
C  ( A SIMPLE BUBBLE SORT )
      IS1=IS+1
      DO 130 J=IS1,IF
      I=J
  120 IF(COUNT(I-1).LE.COUNT(I))GOTO 130
      AV=COUNT(I-1)
      COUNT(I-1)=COUNT(I)
      COUNT(I)=AV
      I=I-1
      IF(I.GT.IS)GOTO  120
  130 CONTINUE
      LA=LA-2
      GOTO 260
C             *******  QUICKSORT  ********
C  SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS
C  THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S
C  HIGHEST ADDRESS.
  140 IY=(IS+IF)/2
      X=COUNT(IY)
      COUNT(IY)=COUNT(IF)
C  THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END
C  OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE
C  OF X .
      K=1
      IFK=IF
C  WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE
C  INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AS NECESSARY,
C  UNTIL THEY MEET .
      DO 160 I=IS,IF
      IF(X.GT.COUNT(I))GOTO 160
      IF(I.GE.IFK)GOTO 170
      COUNT(IFK)=COUNT(I)
      K1=K
      DO 150 K=K1,IFKA
      IFK=IF-K
      IF(COUNT(IFK).GE.X)GOTO 150
      IF(I.GE.IFK)GOTO 180
      COUNT(I)=COUNT(IFK)
      GO TO 160
  150 CONTINUE
      GOTO 170
  160 CONTINUE
C  RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER
C  WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO
C  2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN ANY ELEMENT
C  IN THE SECOND PART, AND THEY MAY NOW BE SORTED INDEPENDENTLY.
  170 COUNT(IFK)=X
      IP=IFK
      GOTO 190
  180 COUNT(I)=X
      IP=I
C  STORE THE LONGER SUBDIVISION IN WORKSPACE.
  190 IF((IP-IS).GT.(IF-IP))GOTO 200
      MARK(LA)=IF
      MARK(LA-1)=IP+1
      IF=IP-1
      GOTO 210
  200 MARK(LA)=IP-1
      MARK(LA-1)=IS
      IS=IP+1
C  FIND THE LENGTH OF THE SHORTER SUBDIVISION.
  210 LNGTH=IF-IS
      IF(LNGTH)230,260,220
C  IF IT CONTAINS MORE THAN ONE ELEMENT STORE IT IN WORKSPACE .
  220 LA=LA+2
      MARK(LA)=IF
      MARK(LA-1)=IS
      GOTO 270
C  IF IT CONTAINS NO ELEMENTS RESELECT THE OTHER SUBDIVISION
C  AND FIND A DIFFERENT TEST NUMBER. NUMBERS WHICH ARE FOUND TO
C  EQUAL THE TEST NUMBER ARE SORTED OUT.
  230 IS=MARK(LA-1)
      IF=MARK(LA)
      IFEND=IF-1
      IK=IS
      DO 240 I=IK,IFEND
      IF(COUNT(I).NE.X)GOTO 250
      IS=IS+1
  240 CONTINUE
      LA=LA-2
      GO TO 260
  250 AV=COUNT(I)
      IY=(IF+IS)/2
      COUNT(I)=COUNT(IY)
      COUNT(IY)=AV
      GOTO 270
  260 IF(LA.LE.0)GOTO 280
C  OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT
      IF=MARK(LA)
      IS=MARK(LA-1)
  270 CONTINUE
  280 RETURN
      END
C######DATE   12 DEC 1990,     COPYRIGHT UKAEA, HARWELL.
C######ALIAS NB01A
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS NB01A
      SUBROUTINE NB01A (K,AZ,BZ,E2,X,Y,MAXIT)
      implicit none
c
      integer k,maxit,j1,it,m,h,j2,j3
      real az,bz,e2,x,y,a,b,ya,x1,y1,x2,yb,u,ytest
c
      SAVE
      IF(K .GT. 0)GO TO 30
C
C     CALCULATE Y(X) AT X=AZ.
      A = AZ
      B = BZ
      X = A
      J1 = 1
      IT = 1
      M = IABS(MAXIT)
   10 K = J1
   20 RETURN
C
C     PRINT X AND Y(X) WHEN REQUESTED.
   30 IF(MAXIT .LE. 0)WRITE(6,40)X,Y
   40 FORMAT(2X,8HNB01A   ,2HX=,E16.7,2X,5HY(X)=,E16.7)
C
C     TEST WHETHER Y(X) IS SUFFICIENTLY SMALL.
      IF( ABS(Y) .GT. E2)GO TO 50
   45 K = 2
      GO TO 20
C
C     BRANCH DEPENDING ON THE VALUE OF J1.
   50 GO TO (60,70,100,170),J1
C
C     CALCULATE Y(X) AT X=BZ.
   60 YA = Y
      X = B
      J1 = 2
      GO TO 20
C
C     TEST WHETHER THE SIGNS OF Y(AZ) AND Y(BZ) ARE DIFFERENT.
   70 IF(YA*Y .LT. 0.)GO TO 120
C
C     BEGIN THE BINARY SUBDIVISION TO SEARCH FOR A BRACKET.
      X1 = A
      Y1 = YA
      J1 = 3
      H = B-A
      J2 = 1
   80 X2 = A+0.5*H
      J3 = 1
C
C     CHECK WHETHER MAXIT FUNCTION VALUES HAVE BEEN CALCULATED.
   90 IT = IT+1
      IF(IT .GE. M)GO TO 10
      X = X2
      GO TO 20
C
C     TEST WHETHER A BRACKET HAS BEEN FOUND.
  100 IF(YA*Y .LT. 0.)GO TO 120
C
C     CONTINUE THE SEARCH FOR A BRACKET.
      IF(J3 .GE. J2)GO TO 110
      A = X
      YA = Y
      X2 = X+H
      J3 = J3+1
      GO TO 90
  110 A = X1
      YA = Y1
      H = 0.5*H
      J2 = J2+J2
      GO TO 80
C
C     AT THIS POINT THE FIRST BRACKET HAS BEEN FOUND.
  120 B = X
      YB = Y
      J1 = 4
C
C     CALCULATE THE NEXT X BY THE SECANT METHOD BASED ON THE BRACKET.
  130 IF( ABS(YA) .LE. ABS(YB))GO TO 140
      X1 = A
      Y1 = YA
      X = B
      Y = YB
      GO TO 150
  140 X1 = B
      Y1 = YB
      X = A
      Y = YA
C
C     USE THE SECANT METHOD BASED ON THE FUNCTION VALUES Y1 AND Y.
  150 U = Y*(X-X1)/(Y-Y1)
  155 X2 = X-U
      IF(X2.EQ.X)GO TO 195
      X1 = X
      Y1 = Y
      YTEST = 0.5*AMIN1( ABS(YA), ABS(YB))
C
C     CHECK THAT X2 IS INSIDE THE INTERVAL (A,B).
      IF((X2-A)*(X2-B) .LT. 0. )GO TO 90
C
C     CALCULATE THE NEXT VALUE OF X BY BISECTION.
  160 X2 = 0.5*(A+B)
      YTEST = 0.
C
C     CHECK WHETHER THE MAXIMUM ACCURACY HAS BEEN ACHIEVED.
      IF((X2-A)*(X2-B))90,45,45
C
C     REVISE THE BRACKET (A,B).
  170 IF(YA*Y .GE. 0.)GO TO 180
      B = X
      YB = Y
      GO TO 190
  180 A = X
      YA = Y
C
C     USE YTEST TO DECIDE THE METHOD FOR THE NEXT VALUE OF X.
  190 IF(YTEST .LE. 0.)GO TO 130
      IF( ABS(Y) -YTEST)150,150,160
  195 IF(U.EQ.0.)GO TO 45
      U=U+U
      GO TO 155
      END
C######DATE   12 DEC 1990,     COPYRIGHT UKAEA, HARWELL.
C######ALIAS GA15AD
C######CALLS FD05 KB05 NB01
C######DATE   24 JUN 1985     COPYRIGHT UKAEA, HARWELL.
C######ALIAS GA15AD
C###### CALLS KB05  NB01
      SUBROUTINE GA15AD(N, KIND, W, IW, IND, LND, X, Y, TT, XD, YD, LP)
      implicit none
c
      integer n,iw,ind,lnd,kind,lp
c
      DOUBLE PRECISION A, B, BIG, C, CUBIC, D, DM, DP, DY, DYN,
     * EIGHT, EPS, FD05AD
      DOUBLE PRECISION ERR, F, FOUR, ONE, R, ROOT, ROOT1, RT, T, THREE
      DOUBLE PRECISION TL, TN, TT, TWO, T1, T2, W, W2, X, XD, XL
      DOUBLE PRECISION XMAX, XMIN, XN, Y, YD, YL, YN, YP, ZERO
      DIMENSION W(IW), X(N), Y(N), XD(N), YD(N), TT(N)
      DIMENSION IND(2,LND)
      DATA ZERO /0D0/, ONE /1D0/,  TWO /2D0/, THREE /3D0/,
     * FOUR /4D0/,  EIGHT /8D0/
c
c     Local variables
c 
      integer i,i1,m,i2,i3,i4,n1,nstat,ir,nrep,nnd,ncrit
      integer ist,inc,k,j,in,in1,il
c
      CUBIC(R) = D + R*(C+R*(B+R*A))
      BIG=FD05AD(5)*.1D0
      EPS=FD05AD(1)*100.D0
C  CHECK N GREATER THAN 4
      IF (N.LT.4) GO TO 430
C  CHECK DATA IS PERIODIC
      IF (X(1).NE.X(N) .OR. Y(1).NE.Y(N)) GO TO 420
C  CHECK KIND=1 OR KIND=2
      IF (KIND.NE.1 .AND. KIND.NE.2) GO TO 440
      IF (KIND.EQ.1) GO TO 20
C  CHECK VALUES OF TT IN ORDER (CUBICS ONLY)
      DO 10 I=2,N
        IF (TT(I).LE.TT(I-1)) GO TO 450
   10 CONTINUE
C  CHECK XD AND YD ARE PERIODIC  (CUBICS ONLY)
      IF (XD(1).NE.XD(N) .OR. YD(1).NE.YD(N)) GO TO 460
C
C     PARTITION THE STORAGE IN W
C     IN THE PARAMETRIC CUBIC CASE W(I),W(N+I),I=2,N HOLD THE MINIMUM
C  AND MAXIMUM VALUES OF X IN THE INTERVAL BETWEEN KNOTS I-1 AND I.
   20 I1 = 0
      IF (KIND.NE.1) I1 = 2*N
C     W(I1+I),I=1,2,... HOLD THE STATIONARY VALUES OF Y IN ASCENDING
C   ORDER.
      M = (IW-I1)/4
      I2 = I1 + M
C     W(I2+I),I=1,2,... HOLD THE VALUES OF CAPITAL T AT CRITICAL
C  POINTS IN SEQUENCE ROUND THE CURVE.
      I3 = I2 + M
      I4 = I3 + M
C     W(I3+I),W(I4+I),I=1,2,... HOLD THE MAXIMUM AND MINIMUM X VALUES
C  IN THE BOUNDARY SEGMENT ENDING AT THE ITH CRITICAL POINT.
      N1 = N - 1
      NSTAT = 0
      IF (KIND.NE.1) GO TO 50
C
C     FIND MAXIMA AND MINIMA OF Y IN POLYGONAL CASE
C     THE MAXIMA AND MINIMA OF Y AND CORRESPONDING CAPITAL T VALUES
C  ARE STORED TEMPORARILY IN W(I3+I),W(I4+I),I=1,2,...
      DY = Y(N) - Y(N-1)
      DO 40 I=2,N
        DYN = Y(I) - Y(I-1)
        IF (DY*DYN.GT.ZERO) GO TO 30
        NSTAT = NSTAT + 1
        IF (I3+NSTAT.GT.I4) GO TO 400
        W(I3+NSTAT) = Y(I-1)
        W(I4+NSTAT) = I - 1
   30   DY = DYN
   40 CONTINUE
      GO TO 140
C
C     FIND MAXIMA AND MINIMA OF X AND Y.
C     THE MAXIMA AND MINIMA OF Y AND CORRESPONDING CAPITAL T VALUES
C  ARE STORED TEMPORARILY IN W(I3+I),W(I4+I),I=1,2,...
   50 DO 130 I=2,N
        T1 = I - 1
        T2 = I
        CALL GA15CD(T1, T2, X, TT, XD, XMAX, XMIN, N)
        W(N+I) = XMAX
        W(I) = XMIN
        F = (Y(I)-Y(I-1))/FOUR
        A = (TT(I)-TT(I-1))/EIGHT
        DM = YD(I-1)*A
        DP = YD(I)*A
        A = THREE*(DM+DP-F)
        B = DP - DM
        C = THREE*F - DM - DP
        ROOT = -TWO
        ROOT1 = TWO
        IF (B*B.LE.A*C) GO TO 60
        D = B + DSIGN(DSQRT(B*B-A*C),B)
        ROOT = -C/D
        IF (A.NE.ZERO) ROOT1 = -D/A
C     TEST FOR STATIONARY POINTS AT ENDS OF RANGE.
   60   D = (DABS(A)+TWO*DABS(B)+DABS(C))*EPS
        IF (FOUR*DABS(DP).GT.D) GO TO 80
        IF (DABS(ROOT-ONE).LT.DABS(ROOT1-ONE)) GO TO 70
        ROOT1 = ONE
        GO TO 80
   70   ROOT = ONE
   80   IF (FOUR*DABS(DM).GT.D) GO TO 100
        IF (DABS(ROOT+ONE).LT.DABS(ROOT1+ONE)) GO TO 90
        ROOT1 = -ONE
        GO TO 100
   90   ROOT = -ONE
C     STORE STATIONARY POINTS
  100   A = DM + DP - F
        D = Y(I) - A - B - C
        RT = DMIN1(ROOT,ROOT1)
        DO 120 IR=1,2
          IF (DABS(RT).GT.ONE) GO TO 110
          NSTAT = NSTAT + 1
          IF (NSTAT+I3.GT.I4) GO TO 400
          W(NSTAT+I4) = I + (RT-ONE)/TWO
          W(NSTAT+I3) = CUBIC(RT)
  110     RT = DMAX1(ROOT,ROOT1)
  120   CONTINUE
  130 CONTINUE
C
C     SORT STATIONARY POINTS INTO ASCENDING ORDER.
  140 IND(2,1) = NSTAT
      IF (NSTAT.GT.LND) GO TO 410
      DO 150 I=1,NSTAT
        IND(2,I) = -1
        W(I+I1) = W(I+I3)
  150 CONTINUE
      CALL KB05AD(W(I1+1), NSTAT)
C
C     REMOVE ANY REPEATED STATIONARY POINTS.
      NREP = 0
      DO 160 I=2,NSTAT
        IF (W(I1+I-1).EQ.W(I1+I)) NREP = NREP + 1
        W(I1+I-NREP) = W(I1+I)
  160 CONTINUE
C
C     FIND CRITICAL POINTS AND INDEX THE BOUNDARY INTERVALS.
C     FOR EACH INTERVAL BETWEEN W(I1+I-1) AND W(I1+I),I=2,3,...
C  AN ASSOCIATED LIST OF INTEGERS IS HELD IN IND STARTING AT IND(1,I).
C  THE INTEGERS THEMSELVES ARE IN IND(1,.) AND WITH EACH IS
C  ASSOCIATED THE ADDRESS, HELD IN IND(2,.) OF THE NEXT INTEGER OR
C  ZERO TO INDICATE THE END.  EACH INTEGER K IN THIS LIST POINTS
C  TO AN INTERVAL W(I2+IABS(K)-1) TO W(I2+IABS(K)), WHOSE Y VALUES SPAN
C  THE INTERVAL W(I1+I-1) TO W(I1+I).  Y IS INCREASING OR DECREASING
C  IN THIS INTERVAL ACCORDING AS K IS POSITIVE OR NEGATIVE.
      YL = W(I3+NSTAT)
      TL = W(I4+NSTAT)
      NND = NSTAT - NREP
      NCRIT = 0
      IST = 0
  170 IST = IST + 1
      IF (W(I1+IST).NE.YL) GO TO 170
      DO 310 I=1,NSTAT
        YN = W(I3+I)
        TN = W(I4+I)
        IF (TL.GT.TN) TN = TN + N1
        INC = 1
        IF (YN.LT.YL) INC = -1
  180   NCRIT = NCRIT + 1
        IF (I2+NCRIT.GT.I3) GO TO 400
        W(I2+NCRIT) = W(I4+I)
        IF (YN.EQ.YL) GO TO 300
        IST = IST + INC
        YP = W(I1+IST)
        IF (YP.EQ.YN) GO TO 250
        K = TL
  190   J = K
        K = K + 1
        IF (J.GE.N) J = J - N1
        IF ((Y(J+1)-YP)*INC.LE.ZERO .AND. K.LT.TN) GO TO 190
        T1 = DMAX1(J+ZERO,TL)
        T2 = DMIN1(J+ONE,TN)
        IF (KIND.EQ.1) GO TO 240
C     FIND THE INTERCEPT OF Y=YP WITH THE  SPLINE IN Y
C    IN THE INTERVAL T1,T2 (ASSUMED TO LIE IN ONE KNOT INTERVAL).
        J = T1
        A = (TT(J+1)-TT(J))/EIGHT
        DM = YD(J)*A
        DP = YD(J+1)*A
        F = (Y(J+1)-Y(J))/FOUR
        A = DM + DP - F
        B = DP - DM
        C = THREE*F - DM - DP
        D = Y(J+1) - A - B - C
        K = 0
        T = DMAX1(DABS(T1),DABS(T2))*TWO
        ERR = EPS*(DABS(D)+DABS(C)*(ONE+T)+DABS(B)*(ONE+T+T)+DABS(A)*
     *   (ONE+T+T+T))
  200   CALL NB01AD(K, T1, T2, ERR, T, R, 100)
        GO TO (210, 230, 220, 220), K
  210   R = CUBIC((T-J)*TWO-ONE) - YP
        GO TO 200
  220   WRITE (LP,99999)
        IND(1,1) = 4
        GO TO 390
  230   W(I2+NCRIT) = T
        GO TO 250
  240   W(I2+NCRIT) = J + 1 + (YP-Y(J+1))/(Y(J+1)-Y(J))
C     INSERT CRITICAL INTERVAL IN APPROPRIATE LIST.
  250   J = IST - (INC-1)/2
  260   IF (IND(2,J)) 290, 280, 270
  270   J = IND(2,J)
        GO TO 260
  280   IND(2,J) = NND + 1
        NND = NND + 1
        IF (NND.GT.LND) GO TO 410
        J = NND
  290   IND(1,J) = NCRIT*INC
        IND(2,J) = 0
        IF (YP.NE.YN) GO TO 180
  300   TL = W(I4+I)
        YL = YN
  310 CONTINUE
      NSTAT = NSTAT - NREP
C
C     FIND MAXIMUM AND MINIMUM X VALUES FOR EACH INTERVAL BETWEEN
C  CRITICAL POINTS.
      TL = W(I2+NCRIT)
      DO 380 I=1,NCRIT
        TN = W(I2+I)
        IN = TN
        IN1 = IN + 1
        IF (IN1.GT.N) IN1 = 2
        XN = X(IN) + (TN-IN)*(X(IN1)-X(IN))
        IF (TL.GT.TN) TN = TN + N1
        IL = TL
        IN = TN
        IF (KIND.NE.1) GO TO 330
        XL = X(IL) + (TL-IL)*(X(IL+1)-X(IL))
        W(I3+I) = DMAX1(XL,XN)
        W(I4+I) = DMIN1(XL,XN)
        IL = IL + 1
        IF (IL.GT.IN) GO TO 370
        DO 320 K=IL,IN
          J = K
          IF (J.GE.N) J = J - N1
          W(I3+I) = DMAX1(W(I3+I),X(J))
          W(I4+I) = DMIN1(W(I4+I),X(J))
  320   CONTINUE
        GO TO 370
  330   DO 360 J=IL,IN
          T1 = DMAX1(J+ZERO,TL)
          T2 = DMIN1(J+ONE,TN)
          IF (T2.LE.N) GO TO 340
          T1 = T1 - N1
          T2 = T2 - N1
  340     CALL GA15CD(T1, T2, X, TT, XD, XMAX, XMIN, N)
          IF (J.NE.IL) GO TO 350
          W(I3+I) = XMAX
          W(I4+I) = XMIN
  350     W(I3+I) = DMAX1(W(I3+I),XMAX)
          W(I4+I) = DMIN1(W(I4+I),XMIN)
  360   CONTINUE
  370   TL = W(I2+I)
  380 CONTINUE
      IND(1,1) = -NSTAT
      IND(2,1) = NCRIT
      RETURN
C
C     ERROR RETURNS
  390 WRITE (LP,99998)
      RETURN
  400 WRITE (LP,99997)
      IND(1,1) = 1
      GO TO 390
  410 WRITE (LP,99996)
      IND(1,1) = 2
      GO TO 390
  420 WRITE (LP,99995)
      IND(1,1) = 3
      GO TO 390
  430 WRITE (LP,99994)
      IND(1,1) = 5
      GO TO 390
  440 WRITE (LP,99993)
      IND(1,1) = 6
      GO TO 390
  450 WRITE (LP,99992)
      IND(1,1) = 7
      GO TO 390
  460 WRITE (LP,99991)
      IND(1,1) = 8
      GO TO 390
99999 FORMAT (32X, 15HOF NB01 FAILURE)
99998 FORMAT (31H+ERROR RETURN FROM GA15 BECAUSE)
99997 FORMAT (32X, 15HIW IS TOO SMALL)
99996 FORMAT (32X, 16HLND IS TOO SMALL)
99995 FORMAT (32X, 24HTHE DATA IS NOT PERIODIC)
99994 FORMAT (32X, 14HN IS TOO SMALL)
99993 FORMAT (32X, 24HKIND NOT EQUAL TO 1 OR 2)
99992 FORMAT (32X, 14HT NOT IN ORDER)
99991 FORMAT (32X, 24HXD OR YD IS NOT PERIODIC)
      END
      SUBROUTINE GA15BD(XP, YP, RESULT, N, KIND, W, IW, IND, LND, X, Y,
     * TT, XD, YD, LP)
      implicit none
c
      integer n,iw,ind,lnd,lp,kind
c
      DOUBLE PRECISION A, B, BIG, C, CUBIC, D, DIST, DM, DP, EIGHT, EPS
      DOUBLE PRECISION ERR, F, FOUR, ONE, R, RESL, RESULT
      DOUBLE PRECISION T, THREE, TT, TWO, T1, T2, W, X, XD, XDIST
      DOUBLE PRECISION XP, Y, YD, YP, ZERO, FD05AD
      DIMENSION W(IW), X(N), Y(N), XD(N), YD(N), TT(N)
      DIMENSION IND(2,LND)
      DATA ZERO /0D0/, ONE /1D0/,  TWO /2D0/, THREE /3D0/,
     * FOUR /4D0/,  EIGHT /8D0/
c
c     Local variables:
c
      integer i1,m,i2,i3,i4,n1,nstat,ncrit,j1,ju,jt,nl,j,iu,il,it,k,jl
c
      CUBIC(R) = D + R*(C+R*(B+R*A))
      BIG=FD05AD(5)*.1D0
      EPS=FD05AD(1)*100.0D0
C     PARTITION THE STORAGE IN W
      I1 = 0
      IF (KIND.NE.1) I1 = 2*N
      M = (IW-I1)/4
      I2 = I1 + M
      I3 = I2 + M
      I4 = I3 + M
      N1 = N - 1
      NSTAT = -IND(1,1)
      NCRIT = IND(2,1)
C
C     FIND THE INTERVAL BETWEEN STATIONARY POINTS.
C     W(I1+JL).LT.YP.LE.W(I1+JU)
      JL = 0
      JU = NSTAT + 1
   10 IF (JU.EQ.JL+1) GO TO 30
      JT = (JL+JU)/2
      IF (YP.GT.W(I1+JT)) GO TO 20
      JU = JT
      GO TO 10
   20 JL = JT
      GO TO 10
   30 DIST = BIG
      RESULT = -BIG
      IF (JU.EQ.1) GO TO 220
      IF (JU.GT.NSTAT) GO TO 230
C
C     COUNT THE BOUNDARY SEGMENTS TO THE LEFT IN NL.  XDIST IS SET
C  POSITIVE OR NEGATIVE ACCORDING AS THE BOUNDARY SEGMENT IS TO THE LEFT
C  OR RIGHT.  ABS(XDIST) IS A LOWER BOUND ON DISTANCE TO THE BOUNDARY
C  AND IS ZERO FOR POINTS ON THE BOUNDARY.
   40 NL = 0
   50 J = IND(1,JU)
      F = J
      J = IABS(J)
C     TEST AGAINST MAXIMUM AND MINIMUM X VALUES FOR THE BOUNDARY
C  SEGMENT.
      XDIST = XP - W(I3+J)
      IF (XDIST.GT.ZERO) GO TO 200
      XDIST = XP - W(I4+J)
      IF (XDIST.LT.ZERO) GO TO 210
C     PERFORM A BINARY SEARCH TO FIND KNOT INTERVAL INVOLVED.
      IU = W(I2+J) + ONE
      IF (IU.EQ.N+1) IU = N
      IF (J.EQ.1) GO TO 60
      IL = W(I2+J-1)
      GO TO 70
   60 IL = W(I2+NCRIT)
   70 IF (IL.GE.IU) IL = IL - N1
      IF (IU.EQ.IL+1) GO TO 90
      IT = (IL+IU)/2
      IF (IT.LE.1) IT = IT + N1
      IF ((YP-Y(IT))*F.GT.ZERO) GO TO 80
      IU = IT
      IL = IL + N1
      GO TO 70
   80 IL = IT
      GO TO 70
C  CHECK FOR VERTEX
   90 XDIST=ZERO
      IF(XP.EQ.X(IU).AND.YP.EQ.Y(IU)) GO TO 190
      IF (KIND.NE.1) GO TO 100
C     PERFORM DIRECT TEST IN THE POLYGONAL CASE.
      XDIST = XP - (YP-Y(IL))*(X(IU)-X(IL))/(Y(IU)-Y(IL)) - X(IL)
      GO TO 190
C     TEST AGAINST MAXIMUM AND MINIMUM X VALUES IN KNOT INTERVAL.
  100 XDIST = XP - W(N+IU)
      IF (XDIST.GT.ZERO) GO TO 200
      XDIST = XP - W(IU)
      IF (XDIST.LT.ZERO) GO TO 210
C     FIND ACTUAL POINT OF INTERSECTION WITH THE BOUNDARY.
      T = W(I2+NCRIT)
      IF (J.NE.1) T = W(I2+J-1)
      T1 = IL
      IF (T.GT.IL .AND. T.LT.IU) GO TO 110
      XDIST = XP - X(IL)
      IF (YP-Y(IL)) 120, 190, 120
  110 T1 = T
  120 T = W(I2+J)
      T2 = IU
      IF (T.GT.IL .AND. T.LT.IU) GO TO 130
      XDIST = XP - X(IU)
      IF (YP.EQ.Y(IU)) GO TO 190
      GO TO 140
  130 T2 = T
C
C     FIND THE INTERCEPT OF Y=YP WITH THE  SPLINE IN Y
C    IN THE INTERVAL T1,T2 (ASSUMED TO LIE IN ONE KNOT INTERVAL).
  140 J = T1
      A = (TT(J+1)-TT(J))/EIGHT
      DM = YD(J)*A
      DP = YD(J+1)*A
      F = (Y(J+1)-Y(J))/FOUR
      A = DM + DP - F
      B = DP - DM
      C = THREE*F - DM - DP
      D = Y(J+1) - A - B - C
      K = 0
      T = DMAX1(DABS(T1),DABS(T2))*TWO
      ERR = EPS*(DABS(D)+DABS(C)*(ONE+T)+DABS(B)*(ONE+T+T)+DABS(A)*
     * (ONE+T+T+T))
  150 CALL NB01AD(K, T1, T2, ERR, T, R, 100)
      GO TO (160, 180, 170, 170), K
  160 R = CUBIC((T-J)*TWO-ONE) - YP
      GO TO 150
  170 WRITE (LP,99999)
      IND(1,1) = 4
      WRITE (LP,99998)
      RETURN
  180 A = (TT(J+1)-TT(J))/EIGHT
      DM = XD(J)*A
      DP = XD(J+1)*A
      F = (X(J+1)-X(J))/FOUR
      A = DM + DP - F
      B = DP - DM
      C = THREE*F - DM - DP
      D = X(J+1) - A - B - C
      XDIST = XP - CUBIC((T-J)*TWO-ONE)
C     TEST XDIST
  190 IF (XDIST) 210, 240, 200
  200 NL = NL + 1
  210 DIST = DMIN1(DIST,DABS(XDIST))
      JU = IND(2,JU)
      IF (JU.NE.0) GO TO 50
C
C     RETURN OR HANDLE THE SPECIAL CASE WHERE YP EQUALS A STATIONARY Y
C  VALUE.
      RESULT = DIST
      IF ((NL/2)*2.EQ.NL) RESULT = -DIST
      IF (YP.EQ.W(I1+JL)) GO TO 250
      JU = JL + 1
  220 IF (YP.EQ.W(I1+JU)) GO TO 260
  230 RETURN
  240 RESULT = ZERO
      RETURN
  250 IF (RESULT*RESL) 240, 230, 230
  260 RESL = -ONE
      IF (JU.EQ.NSTAT) GO TO 250
      JL = JU
      JU = JL + 1
      IF (RESULT.GT.ZERO) RESL = ONE
      GO TO 40
99999 FORMAT (32X, 15HOF NB01 FAILURE)
99998 FORMAT (31H+ERROR RETURN FROM GA15 BECAUSE)
      END
      SUBROUTINE GA15CD(T1, T2, X, TT, XD, XMAX, XMIN, N)
      implicit none
c
      integer n
c
      DOUBLE PRECISION A, B, C, CUBIC, D, DM, DP, E, EIGHT
      DOUBLE PRECISION F, FOUR, ONE, ROOT, R, R1
      DOUBLE PRECISION R2, THREE, TT, TWO, T1, T2, X, XD
      DOUBLE PRECISION XMAX, XMIN, XX, X1, X2, ZERO
      DIMENSION X(N), XD(N), TT(N)
      DATA ZERO /0D0/, ONE /1D0/, TWO /2D0/, THREE /3D0/, FOUR /4D0/,
     * EIGHT /8D0/
c
c     Local variables:
c
      integer l,ir 
c
      CUBIC(R) = D + R*(C+R*(B+R*A))
C
C     THE FOLLOWING INSTRUCTIONS FIND THE MAXIMUM AND MINIMUM VALUES OF
C  X IN THE INTERVAL T1,T2 (ASSUMED TO LIE IN ONE KNOT INTERVAL).
      L = MIN0(N,IDINT(T1)+1)
      F = (X(L)-X(L-1))/FOUR
      A = (TT(L)-TT(L-1))/EIGHT
      DM = XD(L-1)*A
      DP = XD(L)*A
      A = DM + DP - F
      B = DP - DM
      C = THREE*F - DM - DP
      D = X(L) - A - B - C
      R1 = (T1-L)*TWO + ONE
      R2 = (T2-L)*TWO + ONE
      X1 = CUBIC(R1)
      X2 = CUBIC(R2)
      IF (R1.EQ.-ONE) X1 = X(L-1)
      IF (R2.EQ.ONE) X2 = X(L)
      XMAX = DMAX1(X1,X2)
      XMIN = DMIN1(X1,X2)
      IF (B*B.LE.THREE*A*C) GO TO 30
      E = B + DSIGN(DSQRT(B*B-THREE*A*C),B)
      ROOT = -C/E
      DO 20 IR=1,2
        IF (ROOT.LE.R1 .OR. ROOT.GE.R2) GO TO 10
        XX = CUBIC(ROOT)
        XMAX = DMAX1(XMAX,XX)
        XMIN = DMIN1(XMIN,XX)
   10   IF (A.NE.ZERO) ROOT = -E/A
   20 CONTINUE
   30 RETURN
      END
C######DATE   14 MAY 1991,     COPYRIGHT UKAEA, HARWELL.
C######ALIAS FD05AD
*######DATE   27 FEB 1989     COPYRIGHT UKAEA, HARWELL.
*######ALIAS FD05AD
      DOUBLE PRECISION FUNCTION FD05AD( INUM )
      INTEGER INUM
      DOUBLE PRECISION DC( 5 )
C
C  REAL CONSTANTS (DOUBLE PRECISION ARITHMETIC).
C
C  OBTAINED FROM H.S.L. SUBROUTINE ZE02AM.
C  NICK GOULD AND SID MARLOW, HARWELL, JULY 1988.
C
C  DC(1) THE SMALLEST POSITIVE NUMBER: 1.0 + DC(1) > 1.0.
C  DC(2) THE SMALLEST POSITIVE NUMBER: 1.0 - DC(2) < 1.0.
C  DC(3) THE SMALLEST NONZERO +VE REAL NUMBER.
C  DC(4) THE SMALLEST FULL PRECISION +VE REAL NUMBER.
C  DC(5) THE LARGEST FINITE +VE REAL NUMBER.
C
      DATA DC( 1 ) /    0.222044604925031309D-15 /
      DATA DC( 2 ) /    0.138777878078144569D-16 /
      DATA DC( 3 ) /    0.539760534693402790D-78 /
      DATA DC( 4 ) /    0.539760534693402790D-78 /
      DATA DC( 5 ) /    0.723700557733226210D+76 /
      IF ( INUM .LE. 0 .OR. INUM .GE. 6 ) THEN
         PRINT 2000, INUM
         STOP
      ELSE
         FD05AD = DC( INUM )
      ENDIF
      RETURN
 2000 FORMAT( ' INUM =', I3, ' OUT OF RANGE IN FD05AD.',
     *        ' EXECUTION TERMINATED.' )
      END
C######DATE   12 DEC 1990,     COPYRIGHT UKAEA, HARWELL.
C######ALIAS KB05AD
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS KB05AD
      SUBROUTINE KB05AD(COUNT,N)
      implicit none
c
      integer n
C
C             KB05AD      HANDLES DOUBLE PRECISION VARIABLES
      DOUBLE PRECISION COUNT(1),AV,X
C  THE WORK-SPACE 'MARK' OF LENGTH 50 PERMITS UP TO 2**(50/2) NUMBERS
C  TO BE SORTED. THIS IS MORE THAN THE IBM VIRTUAL MEMORY SPACE
C  WILL HOLD .
c
c     Local variables
c
      INTEGER MARK(50)
      integer m,la,is,if,mloop,ifka,is1,j,i,iy,k,ifk,k1
      integer ip,lngth,ifend,ik
c
C  CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED
      IF(N.EQ.1)GOTO 280
      IF(N.GE.1)GO TO 110
      WRITE(6,100)
  100 FORMAT(///20X,65H ***KB05AD***NO NUMBERS TO BE SORTED ** RETURN TO
     2 CALLING PROGRAM )
      GOTO 280
C  'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER
C  THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.
  110 M=12
C  SET UP INITIAL VALUES.
      LA=2
      IS=1
      IF=N
      DO 270 MLOOP=1,N
C  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE .
      IFKA=IF-IS
      IF((IFKA+1).GT.M)GOTO 140
C********* FINAL SORTING ***
C  ( A SIMPLE BUBBLE SORT )
      IS1=IS+1
      DO 130 J=IS1,IF
      I=J
  120 IF(COUNT(I-1).LE.COUNT(I))GOTO 130
      AV=COUNT(I-1)
      COUNT(I-1)=COUNT(I)
      COUNT(I)=AV
      I=I-1
      IF(I.GT.IS)GOTO  120
  130 CONTINUE
      LA=LA-2
      GOTO 260
C             *******  QUICKSORT  ********
C  SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS
C  THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S
C  HIGHEST ADDRESS.
  140 IY=(IS+IF)/2
      X=COUNT(IY)
      COUNT(IY)=COUNT(IF)
C  THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END
C  OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE
C  OF X .
      K=1
      IFK=IF
C  WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE
C  INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AS NECESSARY,
C  UNTIL THEY MEET .
      DO 160 I=IS,IF
      IF(X.GT.COUNT(I))GOTO 160
      IF(I.GE.IFK)GOTO 170
      COUNT(IFK)=COUNT(I)
      K1=K
      DO 150 K=K1,IFKA
      IFK=IF-K
      IF(COUNT(IFK).GE.X)GOTO 150
      IF(I.GE.IFK)GOTO 180
      COUNT(I)=COUNT(IFK)
      GO TO 160
  150 CONTINUE
      GOTO 170
  160 CONTINUE
C  RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER
C  WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO
C  2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN ANY ELEMENT
C  IN THE SECOND PART, AND THEY MAY NOW BE SORTED INDEPENDENTLY.
  170 COUNT(IFK)=X
      IP=IFK
      GOTO 190
  180 COUNT(I)=X
      IP=I
C  STORE THE LONGER SUBDIVISION IN WORKSPACE.
  190 IF((IP-IS).GT.(IF-IP))GOTO 200
      MARK(LA)=IF
      MARK(LA-1)=IP+1
      IF=IP-1
      GOTO 210
  200 MARK(LA)=IP-1
      MARK(LA-1)=IS
      IS=IP+1
C  FIND THE LENGTH OF THE SHORTER SUBDIVISION.
  210 LNGTH=IF-IS
      IF(LNGTH)230,260,220
C  IF IT CONTAINS MORE THAN ONE ELEMENT STORE IT IN WORKSPACE .
  220 LA=LA+2
      MARK(LA)=IF
      MARK(LA-1)=IS
      GOTO 270
C  IF IT CONTAINS NO ELEMENTS RESELECT THE OTHER SUBDIVISION
C  AND FIND A DIFFERENT TEST NUMBER. NUMBERS WHICH ARE FOUND TO
C  EQUAL THE TEST NUMBER ARE SORTED OUT.
  230 IS=MARK(LA-1)
      IF=MARK(LA)
      IFEND=IF-1
      IK=IS
      DO 240 I=IK,IFEND
      IF(COUNT(I).NE.X)GOTO 250
      IS=IS+1
  240 CONTINUE
      LA=LA-2
      GO TO 260
  250 AV=COUNT(I)
      IY=(IF+IS)/2
      COUNT(I)=COUNT(IY)
      COUNT(IY)=AV
      GOTO 270
C  FIND IF SORTING IS COMPLETED.
  260 IF(LA.LE.0)GOTO 280
C  OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT
      IF=MARK(LA)
      IS=MARK(LA-1)
  270 CONTINUE
  280 RETURN
      END
C######DATE   12 DEC 1990,     COPYRIGHT UKAEA, HARWELL.
C######ALIAS NB01AD
C######DATE   01 JAN 1984     COPYRIGHT UKAEA, HARWELL.
C######ALIAS NB01AD
      SUBROUTINE NB01AD(K,AZ,BZ,E2,X,Y,MAXIT)
      implicit none
      integer k,maxit
      DOUBLE PRECISION A,B,H,U,X,Y,AZ,BZ,E2,X1,X2,YA,YB,Y1,YTEST
c
c     Local variables 
c
      integer j1,it,m,j2,j3
c
      SAVE
c

      IF(K .GT. 0)GO TO 30
C
C     CALCULATE Y(X) AT X=AZ.
      A = AZ
      B = BZ
      X = A
      J1 = 1
      IT = 1
      M = IABS(MAXIT)
   10 K = J1
   20 RETURN
C
C     PRINT X AND Y(X) WHEN REQUESTED.
   30 IF(MAXIT .LE. 0)WRITE(6,40)X,Y
   40 FORMAT(2X,8HNB01AD  ,2HX=,D24.16,2X,5HY(X)=,D24.16)
C
C     TEST WHETHER Y(X) IS SUFFICIENTLY SMALL.
      IF(DABS(Y) .GT. E2)GO TO 50
   45 K = 2
      GO TO 20
C
C     BRANCH DEPENDING ON THE VALUE OF J1.
   50 GO TO (60,70,100,170),J1
C
C     CALCULATE Y(X) AT X=BZ.
   60 YA = Y
      X = B
      J1 = 2
      GO TO 20
C
C     TEST WHETHER THE SIGNS OF Y(AZ) AND Y(BZ) ARE DIFFERENT.
   70 IF(YA*Y .LT.0.0D0)GO TO 120
C
C     BEGIN THE BINARY SUBDIVISION TO SEARCH FOR A BRACKET.
      X1 = A
      Y1 = YA
      J1 = 3
      H = B-A
      J2 = 1
   80 X2 = A+0.5D0*H
      J3 = 1
C
C     CHECK WHETHER MAXIT FUNCTION VALUES HAVE BEEN CALCULATED.
   90 IT = IT+1
      IF(IT .GE. M)GO TO 10
      X = X2
      GO TO 20
C
C     TEST WHETHER A BRACKET HAS BEEN FOUND.
  100 IF(YA*Y .LT.0.0D0)GO TO 120
C
C     CONTINUE THE SEARCH FOR A BRACKET.
      IF(J3 .GE. J2)GO TO 110
      A = X
      YA = Y
      X2 = X+H
      J3 = J3+1
      GO TO 90
  110 A = X1
      YA = Y1
      H = 0.5D0*H
      J2 = J2+J2
      GO TO 80
C
C     AT THIS POINT THE FIRST BRACKET HAS BEEN FOUND.
  120 B = X
      YB = Y
      J1 = 4
C
C     CALCULATE THE NEXT X BY THE SECANT METHOD BASED ON THE BRACKET.
  130 IF(DABS(YA) .LE.DABS(YB))GO TO 140
      X1 = A
      Y1 = YA
      X = B
      Y = YB
      GO TO 150
  140 X1 = B
      Y1 = YB
      X = A
      Y = YA
C
C     USE THE SECANT METHOD BASED ON THE FUNCTION VALUES Y1 AND Y.
  150 U = Y*(X-X1)/(Y-Y1)
  155 X2 = X-U
      IF(X2.EQ.X)GO TO 195
      X1 = X
      Y1 = Y
      YTEST =.5D0*DMIN1(DABS(YA),DABS(YB))
C
C     CHECK THAT X2 IS INSIDE THE INTERVAL (A,B).
      IF((X2-A)*(X2-B) .LT. 0.0D0)GO TO 90
C
C     CALCULATE THE NEXT VALUE OF X BY BISECTION.
  160 X2 = 0.5D0*(A+B)
      YTEST = 0.D0
C
C     CHECK WHETHER THE MAXIMUM ACCURACY HAS BEEN ACHIEVED.
      IF((X2-A)*(X2-B))90,45,45
C
C     REVISE THE BRACKET (A,B).
  170 IF(YA*Y .GE.0.0D0)GO TO 180
      B = X
      YB = Y
      GO TO 190
  180 A = X
      YA = Y
C
C     USE YTEST TO DECIDE THE METHOD FOR THE NEXT VALUE OF X.
  190 IF(YTEST .LE.0.0D0)GO TO 130
      IF(DABS(Y) -YTEST)150,150,160
  195 IF(U.EQ.0.0D0)GO TO 45
      U = U+U
      GO TO 155
      END
