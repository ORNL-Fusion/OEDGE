C
C
      SUBROUTINE FIVEC(ER,B,IFLAG,P)

C  VECTORIZED VERSION OF FUNCTION FI FOR GAUSS MEHLER QUADRATURE
C  EVALUATE EFFECTIVE POTENTIAL FUNCTION AT AR(I),I=1,NFI
C  NOTE: NFI.LE.128 IS NOT CHECKED, BUT USED
C  RETURN FI(AR(I)) IN THE ARRAY AFI(I),I=1,NFI
C     --------------
C  IFLAG=1:  H+ + H
C  IFLAG=2:  H+ + NOBLE GASES , H+ + H2, HE+ + HE
C

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: RI, G1, SS, EX, EX2, RR, RLIM, REFF, EMRR, G2, V,
     .          R, R2, B2, RLOW, RR1, G, RMI, EPS
      INTEGER :: IFI

      REAL(DP) :: AR(128), AFI(128)
      INTEGER :: NFI
      COMMON /CFI/ AR,AFI,NFI

      B2=B*B
C
      IF(IFLAG.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.D-7
        DO 1 IFI=1,NFI
          R=AR(IFI)
          R2=R*R
C  FIND V=V(R)
          IF (R.GT.160.D0) THEN
            V=0.D0
          ELSEIF (R.LT.RLOW) THEN
            R2=RLOW*RLOW
            EX=EXP(-RLOW)
            EX2=EX*EX
            SS=(1.+RLOW+R2/3.)*EX-1.
            RI=1./RLOW
            V=27.211*((RI-(1.+RI)*EX2-(1.+RLOW)*EX)/SS+RI)
          ELSE
            EX=EXP(-R)
            EX2=EX*EX
            SS=(1.+R+R2/3.)*EX-1.
            RI=1./R
            V=27.211*((RI-(1.+RI)*EX2-(1.+R)*EX)/SS+RI)
          ENDIF
          AFI(IFI)=1.-V/ER-B2/R2
1       CONTINUE
C
      ELSEIF(IFLAG.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES, (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1./P(4)
        DO 2 IFI=1,NFI
          R=AR(IFI)
          R2=R*R
          RR=R*RMI
          EMRR=1.-RR
          IF (RR.LT.1) THEN
            G=G1
          ELSE
            G=G1*G2
          ENDIF
          REFF=-G*EMRR
C
          IF (REFF.GT.160.D0) THEN
            V=0.D0
          ELSE
            EX=EXP(-REFF)
            EX2=EX*EX
            V=EPS*(EX2-(EX+EX))
          ENDIF
          AFI(IFI)=1.-V/ER-B2/R2
2       CONTINUE
C
      ELSE
        WRITE (iunout,*) 'ERROR IN FUNCTION FIVEC. IFLAG INVALID.'
        WRITE (iunout,*) 'IFLAG = ',IFLAG
        CALL EXIT_OWN(1)
      ENDIF
C
      RETURN
      END
