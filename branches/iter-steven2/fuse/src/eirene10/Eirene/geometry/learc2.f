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
