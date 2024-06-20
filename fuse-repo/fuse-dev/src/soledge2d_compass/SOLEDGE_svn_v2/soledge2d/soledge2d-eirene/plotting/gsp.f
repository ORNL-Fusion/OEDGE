C
C
      SUBROUTINE EIRENE_GSP (P1,P2,R1,R2,Q1,Q2,S1,S2,XLA,XMU,EPS)
 
      USE EIRMOD_PRECISION
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: P1, P2, R1, R2, Q1, Q2, S1, S2, EPS
      REAL(DP), INTENT(OUT) :: XLA, XMU
      REAL(DP) :: S3, R3, PK1, PK2, PK3, XNO
 
      S3=0.
      R3=0.
      PK1=R2*S3-R3*S2
      PK2=R3*S1-R1*S3
      PK3=R1*S2-R2*S1
      XNO=SQRT(PK1*PK1+PK2*PK2+PK3*PK3)
      IF (XNO.LT.EPS) THEN
       XLA=10.
       XMU=10.
      ELSE IF (ABS(R1).LT.EPS) THEN
       XLA=(P1-Q1)/S1
       XMU=-(P2-Q2-S2*XLA)/R2
      ELSE
       XLA=(P2-Q2)/S2
       XMU=-(P1-Q1-S1*XLA)/R1
      ENDIF
      RETURN
      END
