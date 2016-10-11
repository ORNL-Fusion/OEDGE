C
C
C*DK XCORDF
      SUBROUTINE XCORDF(RQ,E,ELQ,AI,BI,X,XS,Y,YS)
C
C   INTERSECTION OF LINE: AI*Y+BI*X=0 WITH ELLIPSE (X-E)**2+Y*Y/ELQ=RQ
C   (X,Y)=INTERSECTION-POINT IS RETURNED
C   IF NO INTERSECTION, X IS SET TO 1.D30
C   THE SOLUTION WITH LARGER X IS (X,Y), THE OTHER IS (XS,YS)
C
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: RQ, E, ELQ, AI, BI
      REAL(DP), INTENT(OUT) :: X, XS, Y, YS
      REAL(DP) :: T, A, AA, AAA, B

       IF (AI.EQ.0.) GOTO 2
         T=-BI/AI
         A=T*T/ELQ
         AA=A+1.
         AAA=1./AA
         B=RQ*AA-E*E*A
         IF (B.LE.0.) GOTO 3
1        B=SQRT(B)
         X=AAA*(E+B)
         Y=T*X
         XS=AAA*(E-B)
         YS=T*XS
         RETURN
2      X=0.
       XS=0.
       IF (RQ.LE.E*E) GOTO 3
       B=SQRT((RQ-E*E)*ELQ)
       Y=B
       YS=-B
       RETURN
C    NO INTERSECTION
3      X=1.D30
       XS=1.D30
       Y=0.
       YS=0.
       RETURN
       END
