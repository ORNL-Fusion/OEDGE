C
C
C##       MA20A          28/06/72
C NAME MA20A(R)                  CHECK
      SUBROUTINE MA20A(Q,D,A,R,S,IQ,M,N,TOLER)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IQ, M, N
      REAL(DP), INTENT(IN) :: TOLER
      REAL(DP), INTENT(OUT) :: Q(IQ,*),A(*),D(*),R(*)
      INTEGER, INTENT(OUT) :: S(*) 
      REAL(DP) :: SUM, B, PIVOT, BIG
      REAL(DP) :: MIN,MAX
      INTEGER :: OUT, KL, KR, IN, L, K, KOUNT, M2, I, J, N1,
     .           N2, M1
      LOGICAL :: STAGE, TEST
C  ***BIG MUST BE SET EQUAL TO ANY VERY LARGE REAL CONSTANT.
C  ***ITS VALUE HERE IS APPROPRIATE FOR THE IBM 370.
      DATA BIG /1.D75/
C  ***INITIALIZATION
      M2=M+2
      N2=N+2
      M1=M+1
      N1=N+1
      DO 1 J=1,N
      Q(M2,J)=J
    1 A(J)=0.
      DO 3 I=1,M
      Q(I,N2)=N+I
      D(I)=0.
      IF(Q(I,N1).GE.0) GO TO 3
      DO 2 J=1,N2
    2 Q(I,J)=-Q(I,J)
    3 CONTINUE
C  ***COMPUTE MARGINAL COSTS
      DO 5 J=1,N1
      SUM=0.
      DO 4 I=1,M
    4 SUM=SUM+Q(I,J)
    5 Q(M1,J)=SUM
C  ***STAGE I
C  ***DETERMINE VECTOR TO ENTER THE BASIS
      STAGE=.TRUE.
      KOUNT=0
      KR=1
      KL=1
    6 MAX=-1.
      DO 7 J=KR,N
      IF(ABS(Q(M2,J)).GT.N) GO TO 7
      B=ABS(Q(M1,J))
      IF(B.LE.MAX) GO TO 7
      MAX=B
      IN=J
    7 CONTINUE
      IF(Q(M1,IN).GE.0) GO TO 9
      DO 8 I=1,M2
    8 Q(I,IN)=-Q(I,IN)
C  ***DETERMINE VECTOR TO LEAVE THE BASIS
    9 K=0
      DO 10 I=KL,M
      B=Q(I,IN)
      IF(B.LE.TOLER) GO TO 10
      K=K+1
      R(K)=Q(I,N1)/B
      S(K)=I
      TEST=.TRUE.
   10 CONTINUE
   11 IF(K.GT.0) GO TO 12
      TEST=.FALSE.
      GO TO 14
   12 MIN=BIG
      DO 13 I=1,K
      IF(R(I).GE.MIN) GO TO 13
      J=I
      MIN=R(I)
      OUT=S(I)
   13 CONTINUE
      R(J)=R(K)
      S(J)=S(K)
      K=K-1
C  ***CHECK FOR LINEAR DEPENDENCE IN STAGE I
   14 IF(TEST.OR..NOT.STAGE) GO TO 16
      DO 15 I=1,M2
      B=Q(I,KR)
      Q(I,KR)=Q(I,IN)
   15 Q(I,IN)=B
      KR=KR+1
      GO TO 25
   16 IF(TEST) GO TO 17
      Q(M2,N1)=2.
      GO TO 34
   17 PIVOT=Q(OUT,IN)
      IF(Q(M1,IN)-PIVOT-PIVOT.LE.TOLER) GO TO 19
      DO 18 J=KR,N1
      B=Q(OUT,J)
      Q(M1,J)=Q(M1,J)-B-B
   18 Q(OUT,J)=-B
      Q(OUT,N2)=-Q(OUT,N2)
      GO TO 11
C  ***PIVOT ON Q(OUT,IN)
   19 DO 20 J=KR,N1
      IF(J.EQ.IN) GO TO 20
      Q(OUT,J)=Q(OUT,J)/PIVOT
   20 CONTINUE
      DO 22 I=1,M1
      IF(I.EQ.OUT) GO TO 22
      B=Q(I,IN)
      DO 21 J=KR,N1
      IF(J.EQ.IN) GO TO 21
      Q(I,J)=Q(I,J)-B*Q(OUT,J)
   21 CONTINUE
   22 CONTINUE
      DO 23 I=1,M1
      IF(I.EQ.OUT) GO TO 23
      Q(I,IN)=-Q(I,IN)/PIVOT
   23 CONTINUE
      Q(OUT,IN)=1./PIVOT
      B=Q(OUT,N2)
      Q(OUT,N2)=Q(M2,IN)
      Q(M2,IN)=B
      KOUNT=KOUNT+1
      IF(.NOT.STAGE) GO TO 26
C  ***INTERCHANGE ROWS IN STAGE I
      KL=KL+1
      DO 24 J=KR,N2
      B=Q(OUT,J)
      Q(OUT,J)=Q(KOUNT,J)
   24 Q(KOUNT,J)=B
   25 IF(KOUNT+KR.NE.N1) GO TO 6
C  ***STAGE II
      STAGE=.FALSE.
C  ***DETERMINE VECTOR TO ENTER THE BASIS
   26 MAX=-BIG
      DO 28 J=KR,N
      B=Q(M1,J)
      IF(B.GE.0) GO TO 27
      IF(B.GT.-2.) GO TO 28
      B=-B-2.
   27 IF(B.LE.MAX) GO TO 28
      MAX=B
      IN=J
   28 CONTINUE
      IF(MAX.LE.TOLER) GO TO 30
      IF(Q(M1,IN).GT.0) GO TO 9
      DO 29 I=1,M2
   29 Q(I,IN)=-Q(I,IN)
      Q(M1,IN)=Q(M1,IN)-2.
      GO TO 9
C  ***PREPARE OUTPUT
   30 L=KL-1
      DO 32 I=1,L
      IF(Q(I,N1).GE.0) GO TO 32
      DO 31 J=KR,N2
   31 Q(I,J)=-Q(I,J)
   32 CONTINUE
      Q(M2,N1)=0.
      IF(KR.NE.1) GO TO 34
      DO 33 J=1,N
      B=ABS(Q(M1,J))
      IF(B.LE.TOLER.OR.2.-B.LE.TOLER) GO TO 34
   33 CONTINUE
      Q(M2,N1)=1.
   34 DO 37 I=1,M
      K=Q(I,N2)
      B=Q(I,N1)
      IF(K.GT.0) GO TO 35
      K=-K
      B=-B
   35 IF(I.GE.KL) GO TO 36
      A(K)=B
      GO TO 37
   36 K=K-N
      D(K)=B
   37 CONTINUE
      Q(M2,N2)=KOUNT
      Q(M1,N2)=N1-KR
      SUM=0.
      DO 38 I=KL,M
   38 SUM=SUM+Q(I,N1)
      Q(M1,N1)=SUM
C     WRITE (6,*) ' A ',(A(I),I=1,M)
      RETURN
      END
