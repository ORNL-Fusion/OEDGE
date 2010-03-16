C
C
      SUBROUTINE ROTADD(A,AI,ILINI,ILEND)
C
C  TRANSFORM CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES
C  CO-ORDINATES OF A POINT IN OLD SYSTEM (X,Y,Z)
C  CO-ORDINATES OF A POINT IN NEW SYSTEM (X',Y',Z')
C
C    X'         X     X          X'                           -1  T
C    Y' =   A * Y  ;  Y  =  AI * Y' ; SUBROUTINE ASSUMES: AI=A  =A
C    Z'         Z     Z          Z'
C
      USE PRECISION
      USE PARMMOD
      USE COMPRT, ONLY: IUNOUT
      USE CADGEO
      USE CTRCEI

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A(3,3),AI(3,3)
      INTEGER,  INTENT(IN) :: ILINI, ILEND
      REAL(DP) :: P1S, P2S, P3S, A3S, A4S, A5S, A6S, A7S, A8S, A9S, S,
     .            DEL, DETER, T, A1S, A2S
      INTEGER  :: I, J
C
C
      DO 100 I=ILINI,ILEND
C
C   CHANGE COEFFICIENTS FOR ALGEBRAIC EQUATION, TRY TO KEEP INVARIANTS
          S=A4LM(I)+A5LM(I)+A6LM(I)
          DEL=DETER(A4LM(I),A7LM(I)/2._DP,A9LM(I)/2._DP,A7LM(I)/2._DP,
     .              A5LM(I),A8LM(I)/2._DP,A9LM(I)/2._DP,A8LM(I)/2._DP,
     .              A6LM(I))
          T=A5LM(I)*A6LM(I)+A6LM(I)*A4LM(I)+A4LM(I)*A5LM(I)-
     .     (A9LM(I)**2+A8LM(I)**2+A7LM(I)**2)*0.25
C         IF (TRCPLT) THEN
C           WRITE (iunout,*) 'INVARIANTS FROM ROTADD: I= ',I
C           WRITE (iunout,*) 'BEFORE: S,DEL,T= ',S,DEL,T
C           WRITE (iunout,*) ' A0...A9 ',
C    .                  A0LM(I),A1LM(I),A2LM(I),A3LM(I),
C    .                  A4LM(I),A5LM(I),A6LM(I),A7LM(I),A8LM(I),A9LM(I)
C         ENDIF
C
C         A0LM(I)=A0LM(I)
          A1S=A1LM(I)*AI(1,1)+A2LM(I)*AI(2,1)+A3LM(I)*AI(3,1)
          A2S=A1LM(I)*AI(1,2)+A2LM(I)*AI(2,2)+A3LM(I)*AI(3,2)
          A3S=A1LM(I)*AI(1,3)+A2LM(I)*AI(2,3)+A3LM(I)*AI(3,3)
          A1LM(I)=A1S
          A2LM(I)=A2S
          A3LM(I)=A3S
C
          A4S=(A4LM(I)*AI(1,1)+A7LM(I)*AI(2,1))*AI(1,1)+
     .        (A5LM(I)*AI(2,1)+A9LM(I)*AI(3,1))*AI(2,1)+
     .        (A6LM(I)*AI(3,1)+A8LM(I)*AI(1,1))*AI(3,1)
          A5S=(A4LM(I)*AI(1,2)+A7LM(I)*AI(2,2))*AI(1,2)+
     .        (A5LM(I)*AI(2,2)+A9LM(I)*AI(3,2))*AI(2,2)+
     .        (A6LM(I)*AI(3,2)+A8LM(I)*AI(1,2))*AI(3,2)
          A6S=(A4LM(I)*AI(1,3)+A7LM(I)*AI(2,3))*AI(1,3)+
     .        (A5LM(I)*AI(2,3)+A9LM(I)*AI(3,3))*AI(2,3)+
     .        (A6LM(I)*AI(3,3)+A8LM(I)*AI(1,3))*AI(3,3)
C         A6S=S-A4S-A5S
C
          A7S=A4LM(I)*AI(1,1)*AI(1,2)*2+
     .        A5LM(I)*AI(2,1)*AI(2,2)*2+
     .        A6LM(I)*AI(3,1)*AI(3,2)*2+
     .        A7LM(I)*(AI(1,1)*AI(2,2)+AI(1,2)*AI(2,1))+
     .        A8LM(I)*(AI(1,1)*AI(3,2)+AI(1,2)*AI(3,1))+
     .        A9LM(I)*(AI(2,1)*AI(3,2)+AI(2,2)*AI(3,1))
          A8S=A4LM(I)*AI(1,1)*AI(1,3)*2+
     .        A5LM(I)*AI(2,1)*AI(2,3)*2+
     .        A6LM(I)*AI(3,1)*AI(3,3)*2+
     .        A7LM(I)*(AI(1,1)*AI(2,3)+AI(1,3)*AI(2,1))+
     .        A8LM(I)*(AI(1,1)*AI(3,3)+AI(1,3)*AI(3,1))+
     .        A9LM(I)*(AI(2,1)*AI(3,3)+AI(2,3)*AI(3,1))
          A9S=A4LM(I)*AI(1,2)*AI(1,3)*2+
     .        A5LM(I)*AI(2,2)*AI(2,3)*2+
     .        A6LM(I)*AI(3,2)*AI(3,3)*2+
     .        A7LM(I)*(AI(2,1)*AI(2,3)+AI(1,3)*AI(2,2))+
     .        A8LM(I)*(AI(1,2)*AI(3,3)+AI(1,3)*AI(3,2))+
     .        A9LM(I)*(AI(2,2)*AI(3,3)+AI(2,3)*AI(3,2))
          A4LM(I)=A4S
          A5LM(I)=A5S
          A6LM(I)=A6S
          A7LM(I)=A7S
          A8LM(I)=A8S
          A9LM(I)=A9S
C
          S=A4LM(I)+A5LM(I)+A6LM(I)
          DEL=DETER(A4LM(I),A7LM(I)/2._DP,A9LM(I)/2._DP,A7LM(I)/2._DP,
     .              A5LM(I),A8LM(I)/2._DP,A9LM(I)/2._DP,A8LM(I)/2._DP,
     .              A6LM(I))
          T=A5LM(I)*A6LM(I)+A6LM(I)*A4LM(I)+A4LM(I)*A5LM(I)-
     .     (A9LM(I)**2+A8LM(I)**2+A7LM(I)**2)*0.25
C         WRITE (iunout,*) 'AFTER:  S,DEL,T= ',S,DEL,T
C
1         IF (RLB(I).LT.0.) THEN
C   CHANGE COEFFICIENTS IN LINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 10 J=1,ILIN(I)
C             ALIMS(J,I)=ALIMS(J,I)
              A1S=XLIMS(J,I)*AI(1,1)+YLIMS(J,I)*AI(2,1)+
     .            ZLIMS(J,I)*AI(3,1)
              A2S=XLIMS(J,I)*AI(1,2)+YLIMS(J,I)*AI(2,2)+
     .            ZLIMS(J,I)*AI(3,2)
              A3S=XLIMS(J,I)*AI(1,3)+YLIMS(J,I)*AI(2,3)+
     .            ZLIMS(J,I)*AI(3,3)
              XLIMS(J,I)=A1S
              YLIMS(J,I)=A2S
              ZLIMS(J,I)=A3S
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
C             ALIMS0(J,I)=ALIMS0(J,I)
              A1S=XLIMS1(J,I)*AI(1,1)+YLIMS1(J,I)*AI(2,1)+
     .            ZLIMS1(J,I)*AI(3,1)
              A2S=XLIMS1(J,I)*AI(1,2)+YLIMS1(J,I)*AI(2,2)+
     .            ZLIMS1(J,I)*AI(3,2)
              A3S=XLIMS1(J,I)*AI(1,3)+YLIMS1(J,I)*AI(2,3)+
     .            ZLIMS1(J,I)*AI(3,3)
              XLIMS1(J,I)=A1S
              YLIMS1(J,I)=A2S
              ZLIMS1(J,I)=A3S
C
              S=XLIMS2(J,I)+YLIMS2(J,I)+ZLIMS2(J,I)
              A4S=XLIMS2(J,I)*AI(1,1)**2+YLIMS2(J,I)*AI(2,1)**2+
     .            ZLIMS2(J,I)*AI(3,1)**2+
     .            XLIMS3(J,I)*AI(1,1)*AI(2,1)+
     .            YLIMS3(J,I)*AI(1,1)*AI(3,1)+
     .            ZLIMS3(J,I)*AI(2,1)*AI(3,1)
              A5S=XLIMS2(J,I)*AI(1,2)**2+YLIMS2(J,I)*AI(2,2)**2+
     .            ZLIMS2(J,I)*AI(3,2)**2+
     .            XLIMS3(J,I)*AI(1,2)*AI(2,2)+
     .            YLIMS3(J,I)*AI(1,2)*AI(3,2)+
     .            ZLIMS3(J,I)*AI(2,2)*AI(3,2)
              A6S=S-A4S-A5S
C
              A7S=XLIMS2(J,I)*(AI(1,1)*AI(1,2)+AI(1,1)*AI(1,2))+
     .            YLIMS2(J,I)*(AI(2,1)*AI(2,2)+AI(2,1)*AI(2,2))+
     .            ZLIMS2(J,I)*(AI(3,1)*AI(3,2)+AI(3,1)*AI(3,2))+
     .            XLIMS3(J,I)*(AI(1,1)*AI(2,2)+AI(1,2)*AI(2,1))+
     .            YLIMS3(J,I)*(AI(1,1)*AI(3,2)+AI(1,2)*AI(3,1))+
     .            ZLIMS3(J,I)*(AI(2,1)*AI(3,2)+AI(2,2)*AI(3,1))
              A8S=XLIMS2(J,I)*(AI(1,1)*AI(1,3)+AI(1,1)*AI(1,3))+
     .            YLIMS2(J,I)*(AI(2,1)*AI(2,3)+AI(2,1)*AI(2,3))+
     .            ZLIMS2(J,I)*(AI(3,1)*AI(3,3)+AI(3,1)*AI(3,3))+
     .            XLIMS3(J,I)*(AI(1,1)*AI(2,3)+AI(1,3)*AI(2,1))+
     .            YLIMS3(J,I)*(AI(1,1)*AI(3,3)+AI(1,3)*AI(3,1))+
     .            ZLIMS3(J,I)*(AI(2,1)*AI(3,3)+AI(2,3)*AI(3,1))
              A9S=XLIMS2(J,I)*(AI(1,2)*AI(1,3)+AI(1,2)*AI(1,3))+
     .            YLIMS2(J,I)*(AI(2,2)*AI(2,3)+AI(2,2)*AI(2,3))+
     .            ZLIMS2(J,I)*(AI(3,2)*AI(3,3)+AI(3,2)*AI(3,3))+
     .            XLIMS3(J,I)*(AI(2,1)*AI(2,3)+AI(1,3)*AI(2,2))+
     .            YLIMS3(J,I)*(AI(1,2)*AI(3,3)+AI(1,3)*AI(3,2))+
     .            ZLIMS3(J,I)*(AI(2,2)*AI(3,3)+AI(2,3)*AI(3,2))
              XLIMS2(J,I)=A4S
              YLIMS2(J,I)=A5S
              ZLIMS2(J,I)=A6S
              XLIMS3(J,I)=A7S
              YLIMS3(J,I)=A8S
              ZLIMS3(J,I)=A9S
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER: CHANGE TO RLB(I)=-6 OPTION AND DEFINE 6 LINEAR
C                     INEQUALITIES
            ALIMS(1,I)=XLIMS1(1,I)
            XLIMS(1,I)=-1.
            YLIMS(1,I)=0.
            ZLIMS(1,I)=0.
            ALIMS(2,I)=-XLIMS2(1,I)
            XLIMS(2,I)=1.
            YLIMS(2,I)=0.
            ZLIMS(2,I)=0.
            ALIMS(3,I)=YLIMS1(1,I)
            XLIMS(3,I)=0.
            YLIMS(3,I)=-1.
            ZLIMS(3,I)=0.
            ALIMS(4,I)=-YLIMS2(1,I)
            XLIMS(4,I)=0.
            YLIMS(4,I)=1.
            ZLIMS(4,I)=0.
            ALIMS(5,I)=ZLIMS1(1,I)
            XLIMS(5,I)=0.
            YLIMS(5,I)=0.
            ZLIMS(5,I)=-1.
            ALIMS(6,I)=-ZLIMS2(1,I)
            XLIMS(6,I)=0.
            YLIMS(6,I)=0.
            ZLIMS(6,I)=1.
            IF (RLB(I).EQ.1.5) THEN
              WRITE (iunout,*) 'ROTADD: TO BE WRITTEN, RLB=1.5'
            ENDIF
            RLB(I)=-6.
            ILIN(I)=6
            ISCN(I)=0
            GOTO 1
        ELSE
C
C   POINT OPTIONS
          P1S=P1(1,I)*A(1,1)+P1(2,I)*A(1,2)+P1(3,I)*A(1,3)
          P2S=P1(1,I)*A(2,1)+P1(2,I)*A(2,2)+P1(3,I)*A(2,3)
          P3S=P1(1,I)*A(3,1)+P1(2,I)*A(3,2)+P1(3,I)*A(3,3)
          P1(1,I)=P1S
          P1(2,I)=P2S
          P1(3,I)=P3S
          P1S=P2(1,I)*A(1,1)+P2(2,I)*A(1,2)+P2(3,I)*A(1,3)
          P2S=P2(1,I)*A(2,1)+P2(2,I)*A(2,2)+P2(3,I)*A(2,3)
          P3S=P2(1,I)*A(3,1)+P2(2,I)*A(3,2)+P2(3,I)*A(3,3)
          P2(1,I)=P1S
          P2(2,I)=P2S
          P2(3,I)=P3S
          IF (RLB(I).GE.3.) THEN
            P1S=P3(1,I)*A(1,1)+P3(2,I)*A(1,2)+P3(3,I)*A(1,3)
            P2S=P3(1,I)*A(2,1)+P3(2,I)*A(2,2)+P3(3,I)*A(2,3)
            P3S=P3(1,I)*A(3,1)+P3(2,I)*A(3,2)+P3(3,I)*A(3,3)
            P3(1,I)=P1S
            P3(2,I)=P2S
            P3(3,I)=P3S
          ENDIF
          IF (RLB(I).GE.4.) THEN
            P1S=P4(1,I)*A(1,1)+P4(2,I)*A(1,2)+P4(3,I)*A(1,3)
            P2S=P4(1,I)*A(2,1)+P4(2,I)*A(2,2)+P4(3,I)*A(2,3)
            P3S=P4(1,I)*A(3,1)+P4(2,I)*A(3,2)+P4(3,I)*A(3,3)
            P4(1,I)=P1S
            P4(2,I)=P2S
            P4(3,I)=P3S
          ENDIF
          IF (RLB(I).GE.5.) THEN
            P1S=P5(1,I)*A(1,1)+P5(2,I)*A(1,2)+P5(3,I)*A(1,3)
            P2S=P5(1,I)*A(2,1)+P5(2,I)*A(2,2)+P5(3,I)*A(2,3)
            P3S=P5(1,I)*A(3,1)+P5(2,I)*A(3,2)+P5(3,I)*A(3,3)
            P5(1,I)=P1S
            P5(2,I)=P2S
            P5(3,I)=P3S
          ENDIF
          IF (RLB(I).GE.6.) THEN
            P1S=P6(1,I)*A(1,1)+P6(2,I)*A(1,2)+P6(3,I)*A(1,3)
            P2S=P6(1,I)*A(2,1)+P6(2,I)*A(2,2)+P6(3,I)*A(2,3)
            P3S=P6(1,I)*A(3,1)+P6(2,I)*A(3,2)+P6(3,I)*A(3,3)
            P6(1,I)=P1S
            P6(2,I)=P2S
            P6(3,I)=P3S
          ENDIF
        ENDIF
C
100   CONTINUE
C
      RETURN
      END
