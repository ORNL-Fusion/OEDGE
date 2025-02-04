C
C
      SUBROUTINE YSHADD (YSH,ILINI,ILEND)
C
C  SHIFT CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES IN Y DIRECTION
C  OLD ORIGIN: YO=O.    (X,Y,Z)  SYSTEM
C  NEW ORIGIN: YN=-YSH   (X,Y',Z) SYSTEM
C
C    Y'=Y+YSH
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: YSH
      INTEGER, INTENT(IN) :: ILINI, ILEND
      REAL(DP) :: YQ
      INTEGER :: I, J
      YQ=YSH*YSH
C
      DO 100 I=ILINI,ILEND
C
C
C   CHANGE COEFFICIENTS FOR ALGEBRAIC EQUATION
          A0LM(I)=A0LM(I)-A2LM(I)*YSH+A5LM(I)*YQ
          A2LM(I)=A2LM(I)-2.*A5LM(I)*YSH
          A1LM(I)=A1LM(I)-A7LM(I)*YSH
          A3LM(I)=A3LM(I)-A9LM(I)*YSH
          IF (RLB(I).LT.0.) THEN
C   CHANGE COEFFICIENTS IN LINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 10 J=1,ILIN(I)
              ALIMS(J,I)=ALIMS(J,I)-YLIMS(J,I)*YSH
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
              ALIMS0(J,I)=ALIMS0(J,I)-YLIMS1(J,I)*YSH+YLIMS2(J,I)*YQ
              YLIMS1(J,I)=YLIMS1(J,I)-2.*YLIMS2(J,I)*YSH
              XLIMS1(J,I)=XLIMS1(J,I)-XLIMS3(J,I)*YSH
              ZLIMS1(J,I)=ZLIMS1(J,I)-ZLIMS3(J,I)*YSH
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER
            YLIMS1(1,I)=YLIMS1(1,I)+YSH
            YLIMS2(1,I)=YLIMS2(1,I)+YSH
C
          ELSE
C
C   POINT OPTIONS
            P1(2,I)=P1(2,I)+YSH
            P2(2,I)=P2(2,I)+YSH
            IF (RLB(I).GE.3.) P3(2,I)=P3(2,I)+YSH
            IF (RLB(I).GE.4.) P4(2,I)=P4(2,I)+YSH
            IF (RLB(I).GE.5.) P5(2,I)=P5(2,I)+YSH
            IF (RLB(I).GE.6.) P6(2,I)=P6(2,I)+YSH
          ENDIF
C
100   CONTINUE
C
      RETURN
      END
