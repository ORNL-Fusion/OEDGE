C
C
      FUNCTION FI(R,ER,B,IFLAG,P,DFI)

C  EVALUATE EFFECTIVE POTENTIAL FUNCTION FI AT R
C  EVALUATE DFI(R)/DR AT R
C  RETURN FI=FI(R), DFI=DFI(R)/DR
C     --------------
C  IFLAG=1:  H+ + H
C  IFLAG=2:  H+ + NOBLE GASES,  H+ + H2,  HE+ + HE
C

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: R, ER, B, P(*)
      REAL(DP), INTENT(OUT) :: DFI
      INTEGER, INTENT(IN) :: IFLAG

      REAL(DP) :: DSS, U, DU, SS, RI, RIQ, G1, REFF, RLIM, RR1, RR, 
     .            EMRR, G2, EX2, B2, V, DV, EX, RLOW, R2, R3, FI, G, 
     .            RMI, EPS

      B2=B*B
C
      IF(IFLAG.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.D-7
        R2=R*R
        R3=R2*R
C  FIND V=V(R) and DV=DV(R)/DR
        IF (R.GT.160.D0) THEN
          V=0.D0
          DV=0.D0
        ELSEIF (R.LT.RLOW) THEN
          R2=RLOW*RLOW
          R3=R2*RLOW
          EX=EXP(-RLOW)
          EX2=EX*EX
          SS=(1.+RLOW+R2/3.)*EX-1.
          RI=1./RLOW
          V=27.211*((RI-(1.+RI)*EX2-(1.+RLOW)*EX)/SS+RI)
          DV=0.D0
        ELSE
          EX=EXP(-R)
          EX2=EX*EX
          SS=(1.+R+R2/3.)*EX-1.
          RI=1./R
          V=27.211*((RI-(1.+RI)*EX2-(1.+R)*EX)/SS+RI)
          RIQ=RI*RI
          DSS=-R/3.*(1.+R)*EX
          U=RI-(1.+RI)*EX2-(1.+R)*EX
          DU=-RIQ+(RIQ+2.*RI+2.)*EX2+R*EX
          DV=27.211*((DU*SS-U*DSS)/(SS*SS)-RIQ)
        ENDIF
C  FIND FI=FI(R) AND DFI=DFI(R)/DR
        FI=1.-V/ER-B2/R2
        DFI=-DV/ER+2.*B2/R3
C
      ELSEIF (IFLAG.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1./P(4)
        R2=R*R
        R3=R2*R
C
        RR=R*RMI
        EMRR=1.-RR
C
C       G2=1.00+(1.0-G2)*MAX(0.D0,-EMRR)/EMRR
C       REFF=-G1*G2*EMRR
        IF (RR.LT.1.0) THEN
          G=G1
        ELSE
          G=G1*G2
        ENDIF
        REFF=-G*EMRR
C
        IF (REFF.GT.160.D0) THEN
          V=0.D0
          DV=0.D0
        ELSE
          EX=EXP(-REFF)
          EX2=EX*EX
          V=EPS*(EX2-(EX+EX))
          DV=-2.*EPS*RMI*G*(EX2-EX)
        ENDIF
        FI=1.-V/ER-B2/R2
        DFI=-DV/ER+2.*B2/R3
C
      ELSE
        WRITE (6,*) 'ERROR IN FUNCTION FI. IFLAG INVALID.'
        WRITE (6,*) 'IFLAG = ',IFLAG
        CALL EXIT_OWN(1)
      ENDIF
C
      RETURN
      END
