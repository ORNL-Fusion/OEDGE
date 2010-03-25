C EIRENE04 COMPILATION
C ===== SOURCE: eirsrt.f

C
      SUBROUTINE EIRSRT(LSTOP,LTIME,DELTAT,FLUXES,
     .                  B2BRM,B2RD,B2Q,B2VP)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: fluxes(*)
      REAL(DP), INTENT(IN) :: DELTAT, B2BRM, B2RD, B2Q, B2VP
      logical, INTENT(IN) :: lstop, ltime
C
      return
      end
C ===== SOURCE: if0prm.f
C
C
      SUBROUTINE IF0PRM(IUNIN)
      USE PRECISION
      USE PARMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIN

      NCOP=1
      NCPV=1

      END
C ===== SOURCE: infcop.f
*DK INFCOP
      SUBROUTINE INFCOP
      IMPLICIT NONE
c slmod begin
      INTEGER, INTENT(IN) :: I1, I2, I3, IPANU
c
c      INTEGER, INTENT(IN) :: I1, I2, I3
c slmod end
      ENTRY IF0COP
      RETURN
      ENTRY IF1COP
      RETURN
      ENTRY IF2COP(I1)
      RETURN
c slmod begin
      ENTRY IF3COP(I1,I2,I3,IPANU)
      CALL OUTUS1(I1,I2,I3,IPANU)
c
c      ENTRY IF3COP(I1,I2,I3)
c slmod end
      RETURN
      ENTRY IF4COP
      RETURN
      END
C ===== SOURCE: mshproj.f
C
C
      SUBROUTINE MSHPROJ(X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                   NDXA,NR1ST,IY)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1(*),Y1(*),X2(*),Y2(*),X3(*),Y3(*),
     .                      X4(*),Y4(*)
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(IN) :: NDXA,NR1ST,IY
      REAL(DP) :: EPS60, D12, D34, D13, D24, DUX, DUY, DVX, DVY, 
     .          PUPV, PVPV
      INTEGER :: IX, IN
      EPS60 = 1.E-60_DP
C
C
      DO 1 IX=1,NDXA
C
C  CALCULATE THE NORM OF THE VECTORS (POINT2-POINT1),....
C
        D12 = SQRT((X2(IX)-X1(IX))*(X2(IX)-X1(IX))+(Y2(IX)-Y1(IX))*
     .        (Y2(IX)-Y1(IX)))+EPS60
        D34 = SQRT((X4(IX)-X3(IX))*(X4(IX)-X3(IX))+(Y4(IX)-Y3(IX))*
     .        (Y4(IX)-Y3(IX)))+EPS60
        D13 = SQRT((X3(IX)-X1(IX))*(X3(IX)-X1(IX))+(Y3(IX)-Y1(IX))*
     .        (Y3(IX)-Y1(IX)))+EPS60
        D24 = SQRT((X4(IX)-X2(IX))*(X4(IX)-X2(IX))+(Y4(IX)-Y2(IX))*
     .        (Y4(IX)-Y2(IX)))+EPS60
C
C  CALCULATE THE BISSECTING VECTORS, BUT NOT NORMALISED YET
C
        DUX = (X2(IX)-X1(IX))/D12 + (X4(IX)-X3(IX))/D34
        DUY = (Y2(IX)-Y1(IX))/D12 + (Y4(IX)-Y3(IX))/D34
        DVX = (X3(IX)-X1(IX))/D13 + (X4(IX)-X2(IX))/D24
        DVY = (Y3(IX)-Y1(IX))/D13 + (Y4(IX)-Y2(IX))/D24
C
C  CALCULATE THE COMPONENTS OF THE TWO UNIT VECTOR (= PROJECTION RATE)
C
        IN=IY+(IX-1)*NR1ST
        PUX(IN) = DUX/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PUY(IN) = DUY/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PVX(IN) = DVX/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
        PVY(IN) = DVY/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
C
C  ORTHOGONORMALIZE, CONSERVE ORIENTATION (E.SCHMIDT)
C
        PUPV=PUX(IN)*PVX(IN)+PUY(IN)*PVY(IN)
        PVX(IN)=PVX(IN)-PUPV*PUX(IN)
        PVY(IN)=PVY(IN)-PUPV*PUY(IN)
        PVPV=SQRT(PVX(IN)*PVX(IN)+PVY(IN)*PVY(IN))+EPS60
        PVX(IN)=PVX(IN)/PVPV
        PVY(IN)=PVY(IN)/PVPV
C
1     CONTINUE
      RETURN
      END
C ===== SOURCE: neutr.f
C
C
*//NEUTR//
C=======================================================================
C          S U B R O U T I N E   N E U T R
C=======================================================================
      SUBROUTINE NEUTR(KARD,NDIMX,NDIMY,NDIMF,DUMMY,LDMX,LDMY,LDMF,
     .                 LDNS,IS)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KARD,NDIMX,NDIMY,NDIMF,LDMX,LDMY,LDMF,
     .                       LDNS,IS
      INTEGER :: ND1,LIM,IX,IY,III,IF
      REAL(DP), INTENT(IN) :: DUMMY(0:LDMX+1,0:LDMY+1,LDMF,LDNS)
C
      ND1 = NDIMX
      LIM = (ND1/5)*5 - 4
      DO  500  IF = 1,NDIMF
        DO  110  IY = 1,NDIMY
          DO  100  IX = 1,LIM,5
  100     WRITE(KARD,910) (DUMMY(IX-1+III,IY,IF,IS),III = 1,5)
          IF( (LIM+4).EQ.ND1 )   GOTO 110
          WRITE(KARD,910) (DUMMY(IX,IY,IF,IS),IX = LIM+5,ND1)
  110   CONTINUE
  500 CONTINUE
      RETURN
  910 FORMAT(5(E16.8))
*//END NEUTR//
      END
C ===== SOURCE: plasm.f


*//PLASM//
C=======================================================================
C          S U B R O U T I N E   P L A S M
C=======================================================================
      SUBROUTINE PLASM(KARD,NDIMX,NDIMY,NDIMF,N,M,NF,DUMMY)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KARD, NDIMX, NDIMY, NDIMF, N, M, NF
      REAL(DP), INTENT(OUT) :: DUMMY(0:N+1,0:M+1,NF)
      INTEGER :: ND1, LIM, IX, IY, IF, III
      ND1 = NDIMX + 2
      LIM = (ND1/5)*5 - 4
      DO    110  IF = 1,NDIMF
      DO    110  IY = 0,NDIMY+1
      DO    100  IX = 1,LIM,5
100     READ(KARD,910) (DUMMY(-1+IX-1+III,IY,IF),III = 1,5)
        IF( (LIM+4).EQ.ND1 )     GOTO 110
        READ(KARD,910) (DUMMY(-1+IX,IY,IF),IX = LIM+5,ND1)
110   CONTINUE
      RETURN
910   FORMAT(5(E16.8))
*//END PLASM//
      END
C ===== SOURCE: statis_cop.f
C
C
      SUBROUTINE STATIS_COP
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1, X2, X3
      INTEGER, INTENT(IN) :: N1, N2, N3, N4, N5
      LOGICAL, INTENT(IN) :: L1, L2
      ENTRY STATS0_COP
      ENTRY STATS1_COP(N1,N2,N3,N4,N5,L1,L2)
      ENTRY STATS2_COP(X1,X2,X3)
      END
C ===== SOURCE: uptcop.f
C
C
c      SUBROUTINE UPTCOP(X,X2,W,IFLAG)
c      USE PRECISION
c      USE PARMMOD
c      USE COMXS
c      IMPLICIT NONE
c      REAL(DP), INTENT(IN) :: X(MSTOR1,MSTOR2,N2ND+N3RD),
c     .                      X2(NSTORV,N2ND+N3RD), W
c      INTEGER, INTENT(IN) :: IFLAG
c       
c      END
