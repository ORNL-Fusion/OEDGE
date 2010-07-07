C*
C*     COLLISIONAL-RADIATIVE MODEL OF
C*
C*     ATOMIC HYDROGEN
C*
C*
C   ASSUME: SLOWLY EVOLVING SPECIES: H,H+
C   ASSUME: QUASI STEADY STATE OF H*(N) WITH H, H+
C   OUTPUT:
C   R0(..)    : TRAIN OF H* TRAVELING WITH H+
C   R1(..)    : TRAIN OF H* TRAVELING WITH H
C   R_EXT(..) : TRAIN OF H* TRAVELING WITH external source, e.g. radiation trap.
C
C
c   C: elec. impact excitation processes
c   F: elec. impact de-excitation processes  (invers to C, detailed balance)
c   A: spontaneous radiative decay
c   S: ionization
c
c   ALPHA(N): H+    ->  H*(N)  THREEBODY recombination from H+
c                              (invers to S: elect. impact ionization)
c   BETA(N) : H+    ->  H*(N)  radiative rec. from H+
c   C(1,N)  : H(1)  ->  H*(N)  excitation from ground state
C   Q_EXT(N): ???   ->  H*(N)  external source
C
c   reduced pop coeff r0,r1  are per electron. hence: times "densel"
c   for pop0,pop1 - arrays of reduced population coefficients
C*
C***********************************************************************
      SUBROUTINE H_COLRAD (TEMP, DENSEL, Q_EXT, POP0, POP1, POP2,
     .                     ALPCR, SCR, SCR_EXT,
     .                     E_ALPCR, E_SCR, E_SCR_EXT,
     .                     E_ALPCR_T, E_SCR_T,E_SCR_EXT_T)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

C--------- ATOMIC PARAMETER ------------------------------------------
      REAL(DP), INTENT(IN) :: TEMP, DENSEL
      REAL(DP), INTENT(IN) :: Q_EXT(40)
      REAL(DP), INTENT(OUT) ::   ALPCR,    SCR,     SCR_EXT
      REAL(DP), INTENT(OUT) :: E_ALPCR,  E_SCR,   E_SCR_EXT
      REAL(DP), INTENT(OUT) :: E_ALPCR_T,E_SCR_T, E_SCR_EXT_T
      REAL(DP), INTENT(OUT) :: POP0(40), POP1(40), POP2(40)

      REAL(DP), SAVE :: A(40,40), E_AT(40), OSC(40,40)
      REAL(DP), SAVE :: A21SAVE, POP_ESC

      REAL(DP) :: C(40,40),S(40),F(40,40)
     &           ,SAHA(40),BETA(40),ALPHA(40)
      REAL(DP) :: R1(40),R0(40),R_EXT(40)
C
cdr   integer, parameter :: lupa=34, lima=40
      integer, parameter :: lupa=34, lima=34
      INTEGER, SAVE :: IFRST=0
      INTEGER :: IP
c
      IF (LIMA.GT.40.OR.LUPA.GT.LIMA) THEN
        WRITE (iunout,*) 'LIMA, LUPA ??? ',LIMA,LUPA
        CALL EXIT_OWN(1)
      ENDIF
C
C ATOM
c  pop_esc   : lyman alpha population escape factor
c  pop_esc= 1: lyman alpha opt. thin
c  pop_esc= 0: lyman alpha opt. thick
c
c  only once and for all !!
c
      IF (IFRST == 0) THEN
        POP_ESC=1.
        CALL EINSTN(OSC,A,E_AT,40,A21SAVE,POP_ESC)
        IFRST = 1
      END IF

c
c
C ATOM
      CALL CLSAHA(TEMP,SAHA)
      CALL RATCOF(TEMP,OSC,SAHA,C,F,S,ALPHA,BETA)

C
C***********************************************************************
C
C ATOMIC HYDROGEN
C
      CALL POPCOF(DENSEL,SAHA,C,F,S,A,ALPHA,BETA,LUPA,LIMA,R0,R1,
     &            R_EXT,Q_EXT)

C TRAINS OF ELECTRONICALLY EXCITED H
      DO IP=2,LIMA
C
CDR COUPLING TO H+ IONS
        POP0(IP)=R0(IP)*DENSEL
CDR COUPLING TO H ATOMS
        POP1(IP)=R1(IP)*DENSEL
CDR COUPLING TO EXTERNAL SOURCE Q
        POP2(IP)=R_EXT(IP)

      END DO
C  EFFECTIVE IONISATION RATES
      CALL IONREC(C,S,SAHA,A,ALPHA,BETA,R0,R1,DENSEL,LUPA,LIMA,
     &            F,R_EXT,Q_EXT,
     &            ALPCR,SCR,SCR_EXT)
C***********************************************************************
C  EFFECTIVE ELECTRON COOLING RATES
      CALL E_IONREC(C,S,SAHA,A,ALPHA,BETA,R0,R1,DENSEL,LUPA,LIMA,
     &              F,R_EXT,Q_EXT,E_AT,
     &              ALPCR,      SCR,    SCR_EXT,
     &              E_ALPCR,  E_SCR,  E_SCR_EXT,
     &              E_ALPCR_T,E_SCR_T,E_SCR_EXT_T)
C

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      RETURN
      END SUBROUTINE H_COLRAD

C**********************************************************************
      SUBROUTINE EINSTN(F,A,E_AT,LIM,A21SAVE,POP_ESC)
C
C     CALCULATION OF OSCILLATOR STRENGTH AND EINSTEIN COEFFICIENT
C     FOR ATOMIC HYDROGEN
C     E_AT are the energy levels. Ground state is E_AT(1)=0.
C
C     L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972).
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(40,40)
      DIMENSION A(40,40)
      DIMENSION E_AT(40)

      UH=13.595
      DO 100 I=1,LIM
        P=I
        E_AT(I)=UH*(1.0-1.0/P**2)
100   CONTINUE

      DO 101 I=1,LIM-1
      DO 102 J=I+1,LIM
      AI=I
      AJ=J
      X=1.-(AI/AJ)**2
      G=(.9935+.2328/AI-.1296/AI**2)
     *-(.6282-.5598/AI+.5299/AI**2)/(AI*X)
C Feb. 2008: TYPO CORRECTED: .3387 --> .3887 
     *+(.3887-1.181/AI+1.470/AI**2)/(AI*X)**2
      IF(I.NE.1.AND.I.NE.2) GO TO 300
      G=1.0785-.2319/X+.02947/X**2
  200 IF(I.NE.1) GO TO 300
      G=1.1330-.4059/X+.07014/X**2
  300 F(I,J)=2.**6/(3.*SQRT(3.)*3.1416)*(AI/AJ)**3/(2.*AI**2)*G/X**3
      A(J,I)=8.03E9*AI**2/AJ**2*(AI**(-2)-AJ**(-2))**2*F(I,J)
  102 CONTINUE
  101 CONTINUE
cdr  nur a(2-->1) rausnehmen
      A21SAVE=A(2,1)
      A(2,1)=A(2,1)*pop_esc
      RETURN
      END

C***********************************************************************
      SUBROUTINE CLSAHA(TEMP,SAHA)
C
C     SAHA-BOLTZMANN COEFFICIENT FOR ATOMIC HYDROGEN
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SAHA(40)

      TE=TEMP*1.1605E4

      DO 101 I=1,40
      P=I
      UION=13.595/TEMP/P**2
c     if (uion.gt.65) then
c       saha(i)=1.e30
c     else
        SAHA(I)=P**2*EXP(UION)/2.414D15/SQRT(TE**3)
c     endif
  101 CONTINUE
      RETURN
      END

C***********************************************************************
      SUBROUTINE RATCOF(TEMP,OSC,SAHA,C,F,S,ALPHA,BETA)
C
C     RATE COEFFICIENT FOR ATOMIC HYDROGEN
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION OSC(40,40),C(40,40),F(40,40),U(40,40)
      DIMENSION SAHA(40),S(40),ALPHA(40),BETA(40),UION(40)

C     INITIALIZATION
      DO 1 I=1,40
      S(I)=0.0
      ALPHA(I)=0.0
      BETA(I)=0.0
      UION(I)=0.0
      DO 1 J=1,40
      C(I,J)=0.0
      F(I,J)=0.0
      U(I,J)=0.0
    1 CONTINUE

      TE=TEMP*1.1605E4
      UH=13.595

      DO 101 I=1,40
      P=I
  101 UION(I)=13.595/TEMP/P**2
      DO 102 I=1,40
      DO 102 J=1,40
  102 U(I,J)=UION(I)-UION(J)

      IF(TE.GT.5.0E3) THEN
        CALL EXCOFF(U,OSC,TEMP,C,F,S,ALPHA)

        DO 105 I=1,40
  105     ALPHA(I)=S(I)*SAHA(I)

      ELSE
        CALL EXCOFF(U,OSC,TEMP,C,F,S,ALPHA)
        S(1)=ALPHA(1)/SAHA(1)
        DO 19 I=2,40
   19     ALPHA(I)=S(I)*SAHA(I)
      END IF
c
      DO 602 I=1,40
      P=I
      XP=UH/TEMP/P**2
      CALL CLBETA(XP,P,XS)
  602 BETA(I)=5.197D-14*(UH/TEMP)**.5/P*XS
c

      RETURN
      END

C***********************************************************************
      SUBROUTINE EXCOFF(U,OSC,TEMP,C,F,S,ALPHA)
C
C     EXCITATION RATE COEFFICIENT
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(40,40),OSC(40,40),C(40,40),F(40,40)
      DIMENSION S(40),ALPHA(40)
C
      DO 1 I=1,40
      S(I)=0.0
      ALPHA(I)=0.0
      DO 1 J=1,40
      C(I,J)=0.0
    1 F(I,J)=0.0

      TE=TEMP*1.1605E4

C*********  1 -> J
      I=1
      P=I
      DO 100 J=2,40
      Q=J
      CALL COF1N(U(I,J),OSC(I,J),TE,F1,I,J)
      F(J,1)=F1
  100 C(1,J)=Q**2/EXP(U(1,J))*F(J,1)  !/P**2, but P=1 here

C*********  2-10 -> J
      DO 110 I=2,10
      P=I
      DO 110 J=I+1,40
      Q=J
      CALL COFVR(U(I,J),OSC(I,J),TEMP,CV,I,J)
      CALL COFJO(U(I,J),OSC(I,J),TE,CJ,I,J)

      GG=((P-2.)/8.)**0.25
      C(I,J)=(1.-GG)*CJ+GG*CV
110   F(J,I)=P**2/Q**2*EXP(U(I,J))*C(I,J)


C*********  I(>11) -> J

      DO 120 I=11,39
      P=I
      DO 120 J=I+1,40
      Q=J
      CALL COFVR(U(I,J),OSC(I,J),TEMP,CV,I,J)
      C(I,J)=CV
  120 F(J,I)=P**2/Q**2*EXP(U(I,J))*C(I,J)

C*********  S  1 ->  ionization
      I=1

      IF(TE.GT.5.0E3) THEN
      CALL COFJS(TE,S1,I)
      S(1)=S1

      ELSE      !  Te  <= 5000
      CALL COFJS2(TE,AL,I)
      ALPHA(1)=AL
      END IF
C
C*********  S  2-10 ->
      DO 210 I=2,10
      P=I

      CALL COFJS(TE,SJ,I)
      CALL COFVS(TEMP,SV,I)

      GGG=((P-2.)/8.)**0.25
  210 S(I)=(1.-GGG)*SJ+GGG*SV

C*********  S  I(>11) ->
      DO 220 I=11,40
      CALL COFVS(TEMP,SV,I)
  220 S(I)=SV
      RETURN
      END

C***********************************************************************
      SUBROUTINE COF1N(U,OSC,TE,F,I,J)
C
C     C(1,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN 1 -> J
C
C
C     K.SAWADA, K.ERIGUCHI, T.FUJIMOTO
C     J. APPL. PHYS. 73, 8122 (1993).
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION GINT
      EXTERNAL GINT

      P=I
      Q=J
      X=1-P**2/Q**2
      Y=U
C
      R=0.45*X
      IF(J.EQ.2) THEN
      RX=-0.95
      ELSE IF(J.EQ.3) THEN
      RX=-0.95
      ELSE IF(J.EQ.4) THEN
      RX=-0.95
      ELSE IF(J.EQ.5) THEN
      RX=-0.95
      ELSE IF(J.EQ.6) THEN
      RX=-0.94
      ELSE IF(J.EQ.7) THEN
      RX=-0.93
      ELSE IF(J.EQ.8) THEN
      RX=-0.92
      ELSE IF(J.EQ.9) THEN
      RX=-0.91
      ELSE
      RX=-0.9
      END IF
C
      Z=R+Y
      A=2*P**2*OSC/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X-0.603/X**2)
      C1=2*P**2/X
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z

      IF(Y.GT.50.0) THEN
      E1Y=EXP(-Y)*GINT(Y)/Y
      E1Z=EXP(-Z)*GINT(Z)/Z
      E2Y=EXP(-Y)*(1.0-GINT(Y))
      E2Z=EXP(-Z)*(1.0-GINT(Z))
      E1=(1/Y+0.5)*E1Y+RX*(1/Z+0.5)*E1Z
      E2=E2Y/Y+RX*E2Z/Z
      F=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)*P**2/Q**2*EXP(Y)
      ELSE

      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)

      E1Y=-E1Y
      E1Z=-E1Z
      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=(1/Y+0.5)*E1Y+RX*(1/Z+0.5)*E1Z
      E2=E2Y/Y+RX*E2Z/Z

      F=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)*P**2/Q**2*EXP(Y)
      END IF
      RETURN
      END

C***********************************************************************
      SUBROUTINE COFJO(U,OSC,TE,C,I,J)
C
C     C(I,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN I -> J
C
C     L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972).
C
C
C
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT REAL*8(A-H,O-Z)
      P=I
      BN=(4.0-18.63/P+36.24/P**2-28.09/P**3)/P
      Q=J
      X=1-P**2/Q**2
      Y=U
      R=1.94/P**1.57*X
      Z=R+Y
      A=2*P**2*OSC/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X+BN/X**2)
      C1=2*P**2/X
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z
      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)
      IF (ICON.NE.0) GOTO 1000
      E1Y=-E1Y
      E1Z=-E1Z
      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=(1/Y+0.5)*E1Y-(1/Z+0.5)*E1Z
      E2=E2Y/Y-E2Z/Z

      C=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)

      RETURN
 1000 WRITE(iunout,*) 'ERROR IN COFJO        ICON = ',ICON
      STOP
      END

C***********************************************************************
      SUBROUTINE COFVR(U,OSC,TEMP,C,I,J)
C
C     C(I,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN I -> J
C
C     L.VRIENS, A.H.M.SMEETS, PHYS. REV. A 22, 940 (1980).
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      UH=13.595
      P=I
      BN=1.4/P*LOG(P)-0.7/P-0.51/P**2+1.16/P**3-0.55/P**4
      Q=J
      S1=Q-P
      X=1-P**2/Q**2
      A=2*P**2*OSC/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X+BN/X**2)
      DELTA=EXP(-B/A)+0.06*S1**2/P**2/Q
      G1=1+TEMP/UH*P**3
      G2=3+11*S1**2/P**2
      G3=6+1.6*S1*Q+0.3/S1**2+0.8*Q**1.5/S1**0.5*(S1-0.6)
      GAMMA=UH*LOG(G1)*G2/G3
      C1=1.6D-7*TEMP**0.5/(TEMP+GAMMA)*EXP(-U)
      C2=0.3*TEMP/UH+DELTA

      C=C1*(A*LOG(C2)+B)

      RETURN
      END

C***********************************************************************
      SUBROUTINE COFJS(TE,S,I)
C
C     S(I); IONIZATION RATE COEFFICIENT FOR ATOMIC HYDROGEN  I ->
C
C
C
C     K.SAWADA, K.ERIGUCHI, T.FUJIMOTO
C     J. APPL. PHYS. 73, 8122 (1993).
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION G(0:2,40)

      G(0,1)=1.1330
      G(1,1)=-0.4059
      G(2,1)=0.07014
      G(0,2)=1.0785
      G(1,2)=-0.2319
      G(2,2)=0.02947
      DO 350 N=3,40
      G(0,N)=0.9935+0.2328/N-0.1296/N**2
      G(1,N)=-0.6282/N+0.5598/N**2-0.5299/N**3
  350 G(2,N)=0.3887/N**2-1.181/N**3+1.470/N**4

      IF (I.EQ.1) THEN

      P=1.0

      Y=1.57770E5/TE
      R=0.45
      Z=R+Y
      RX=-0.59
      A=0.0
      DO 223 K=0,2
  223 A=A+G(K,1)/(K+3)
      A=A*1.9603*P

      B=0.66667*P**2*(5-0.603)
      C1=2*P**2
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z
      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)

      E1Y=-E1Y
      E1Z=-E1Z
      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=1/Y*E1Y+RX*1/Z*E1Z
      E0Y=EXP(-Y)/Y
      E0Z=EXP(-Z)/Z
      EGY=E0Y-2*E1Y+E2Y
      EGZ=E0Z-2*E1Z+E2Z
      E2=EGY+RX*EGZ

      S=1.093D-10*SQRT(TE)*P**2*Y**2*(A*E1+B*E2)

      ELSE
      P=I
      BN=(4.0-18.63/P+36.24/P**2-28.09/P**3)/P
      Y=1.57770E5/TE/P**2
      R=1.94/P**1.57
      Z=R+Y
      A=0.0
      DO 222 K=0,2
  222 A=A+G(K,I)/(K+3)
      A=A*1.9603*P
      B=0.66667*P**2*(5+BN)
      C1=2*P**2
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z
      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)
      E1Y=-E1Y
      E1Z=-E1Z
      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=1/Y*E1Y-1/Z*E1Z
      E0Y=EXP(-Y)/Y
      E0Z=EXP(-Z)/Z
      EGY=E0Y-2*E1Y+E2Y
      EGZ=E0Z-2*E1Z+E2Z
      E2=EGY-EGZ

      S=1.093D-10*SQRT(TE)*P**2*Y**2*(A*E1+B*E2)
      END IF
      RETURN
      END

C***********************************************************************
      SUBROUTINE COFJS2(TE,W,I)
C
C     S(I); IONIZATION RATE COEFFICIENT FOR ATOMIC HYDROGEN  I ->
C     FOR LOW TEMPERATURE
C
C
C     K.SAWADA, K.ERIGUCHI, T.FUJIMOTO
C     J. APPL. PHYS. 73, 8122 (1993).
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION G(0:2,40)
      DOUBLE PRECISION GINT
      EXTERNAL GINT

      G(0,1)= 1.1330
      G(1,1)=-0.4059
      G(2,1)= 0.07014
      G(0,2)= 1.0785
      G(1,2)=-0.2319
      G(2,2)= 0.02947
      DO 350 N=3,40
      G(0,N)=0.9935+0.2328/N-0.1296/N**2
      G(1,N)=-0.6282/N+0.5598/N**2-0.5299/N**3
  350 G(2,N)=0.3887/N**2-1.181/N**3+1.470/N**4

C     IF (I.EQ.1) THEN
      PP=I
      P=1.0
      Y=1.57770E5/TE
      R=0.45
      Z=R+Y
      RX=-0.59

      A=0.0
      DO 223 K=0,2
  223 A=A+G(K,1)/(K+3)
      A=A*1.9603*P

      B=0.66667*P**2*(5-0.603)
      C1=2*P**2
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z

      E1Y=EXP(-Y)*GINT(Y)/Y
      E1Z=EXP(-Z)*GINT(Z)/Z
      E2Y=EXP(-Y)*(1.0-GINT(Y))
      E2Z=EXP(-Z)*(1.0-GINT(Z))
      E1=1/Y*E1Y+RX*1/Z*E1Z
      E0Y=EXP(-Y)/Y
      E0Z=EXP(-Z)/Z
      EGY=E0Y-2*E1Y+E2Y
      EGZ=E0Z-2*E1Z+E2Z
      E2=EGY+RX*EGZ

      W=1.093D-10*SQRT(TE)*P**2*Y**2*(A*E1+B*E2)*EXP(13.595/TE*1.1605E4)
     */2.414D15/SQRT(TE**3)

      RETURN
      END

C***********************************************************************
      SUBROUTINE COFVS(TEMP,S,I)
C
C     S(I); IONIZATION RATE COEFFICIENT FOR ATOMIC HYDROGEN  I ->
C
C
C     L.VRIENS, A.H.M.SMEETS, PHYS. REV. A 22, 940 (1980).
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)

      P=I
      UI=13.595/TEMP/P**2
      UIZ=UI**2.33+4.38*UI**1.72+1.32*UI
      S=9.56D-6/TEMP**1.5*EXP(-UI)/UIZ

      RETURN
      END

cdr
      double precision FUNCTION GINT(xX)
      REAL*8 Xx
      double precision X,GG(8)
cdr
      DATA GG/0.2677737343,8.6347608925,18.0590169730,8.5733287401,
     *        3.9584960228,21.0996530827,25.6329561486,9.5733223454/
cdr
      x=dble(xx)
cdr
      GINT=(GG(1)+GG(2)*X+GG(3)*X**2+GG(4)*X**3+X**4)/(GG(5)+GG(6)
     *      *X+GG(7)*X**2+GG(8)*X**3+X**4)
      RETURN
      END


C***********************************************************************
      SUBROUTINE CLBETA(XP,P,S)
C
C     RADIATIVE RECOMBINATION RATE COEFFICIENT
C
C
C
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 GAUNT3,PP,XPP,A,B,EPSA,EPSR
      COMMON PP,XPP
      EXTERNAL GAUNT3

      II=INT(P)
      PP=P
      XPP=XP

      A=0.0
      B=20.0
cdr   EPSA=1.0D-5
      EPSA=1.0D-4
      EPSR=1.0D-5
      NMIN=15
      NMAX=511

      CALL  AQC8(A,B,GAUNT3,EPSA,EPSR,NMIN,NMAX,S,ERR,N,ICON)
C     IF (ICON.NE.0) GOTO 1000

C     WRITE (iunout,10) II,S,ERR,N,ICON
C  10 FORMAT(1H ,1I3,2(2X,1PD10.3),2X,I4,2X,I5,/)

      RETURN
C1000 WRITE(iunout,*) 'ERROR IN CLBETA       ICON = ',ICON
C     STOP
      END

C***********************************************************************

      REAL*8 FUNCTION GAUNT3(X)
      REAL*8 PP,XPP,U,B,X
      COMMON PP,XPP
      U=X/XPP
      B=PP
      GAUNT3=(1./(U+1.)+0.1728*(U-1.)/B**(2./3.)/(U+1.)**(5./3.)-0.0496*
     *(U**2+4./3.*U+1.)/B**(4./3.)/(U+1.)**(7./3.))*EXP(-X)
      RETURN
      END

C***********************************************************************
      SUBROUTINE POPCOF(DENSEL,SAHA,C,F,S,A,ALPHA,BETA,LUP,LIM,R0,R1,
     &      R_EXT,Q_EXT)
C
C     SOLUTION OF RATE EQUATION FOR ATOMIC HYDROGEN
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8      C(40,40),F(40,40),A(40,40),W(40,40)
     &         ,SAHA(40),S(40),ALPHA(40),BETA(40),R0(40),R1(40)
     &         ,       Q_EXT(40),R_EXT(40)
     &         ,BLAX(40),VW(40),WA(40,40)
      dimension ip(40)

      DO 201 K=2,LUP-1

        DO 202 L=2,K
cdr stoss bevoelkerung von k von unten
  202     W(K,L)=C(L,K)*DENSEL

cc diagonale
cc entvoelkerung durch stoesse nach unten
        SUMF=0.
        DO 301 I=1,K-1
  301     SUMF=SUMF+F(K,I)
cc entvoelkerung durch stoesse nach oben
        SUMC=0.
        DO 302 I=K+1,LIM
  302     SUMC=SUMC+C(K,I)
cc  spontan nach unten
        SUMA=0.
        DO 303 I=1,K-1
  303     SUMA=SUMA+A(K,I)
cdr entvoelkerung von k: stoesse nach unten, nach oben, ionis, spontan
cdr                      nach unten
        W(K,K)=-(DENSEL*(SUMF+SUMC+S(K))+SUMA)

cc diagonale fertig

        DO 203 L=K+1,LUP
cdr bevoelkerung durch: stoesse von oben, spontan von oben
  203     W(K,L)=DENSEL*F(L,K)+A(L,K)

  201 CONTINUE

cdr k loop finished
cdr: jetzt: dito fuer k=lup, d.h. bevoelkerung von oben entfaellt
      DO 211 L=2,LUP-1

  211 W(LUP,L)=C(L,LUP)*DENSEL

      SUMF=0.
      DO 311 I=1,LUP-1
  311 SUMF=SUMF+F(LUP,I)
      SUMC=0.0
      DO 313 I=LUP+1,LIM
  313 SUMC=SUMC+C(LUP,I)
      SUMA=0.
      DO 312 I=1,LUP-1
  312 SUMA=SUMA+A(LUP,I)

      W(LUP,LUP)=-(DENSEL*(SUMF+SUMC+S(LUP))+SUMA)

C  RECHTE SEITEN:

c  vorbereiten fuer recombination
      DO 550 K=2,LUP
        SUMF=0.0
        DO 500 I=LUP+1,LIM
  500     SUMF=SUMF+F(I,K)*SAHA(I)
        SUMAS=0.0
        DO 501 I=LUP+1,LIM
  501     SUMAS=SUMAS+SAHA(I)*A(I,K)
c
c  matrixelemente: 1/s  (densel*rate-coeff. )
c  rechte seiten : cm**3/s, nicht: 1/s, also fuer elektronendichte=1
c                                      (bzw: stosspartnerdichte =1)
c  geht wg. linearitaet.
c  recombination e + H+ --> H*
        W(K,LUP+1)=-(DENSEL*SUMF+SUMAS+(DENSEL*ALPHA(K)+BETA(K)))
c  ionisation e + H --> H*
        W(K,LUP+2)=-C(1,K)
c  external source: Q_EXT
        W(K,LUP+3)=-Q_EXT(K)
  550 CONTINUE

cdr w besetzt fuer w(i,j) i=2,lup,j=2,lup+3
cdr geht gut, solange lup<38
      DO 3001 II=LUP,LUP+2
cdr reduziere w indices um 1: auf wa: i=1,lup-1,j=1,(lup-1)+3
        DO 402 I=1,LUP-1
        DO 402 J=1,LUP-1+3
  402     WA(I,J)=W(I+1,J+1)
        DO 3000 J=1,LUP-1
          BLAX(J)=WA(J,II)
 3000   CONTINUE

        CALL LAX(WA,40,LUP-1,  BLAX,0.0,1,IS,VW,IP,ICON)
c

        DO 3010 J=1,LUP-1
          IF(II.EQ.LUP) THEN
coupling to H+
            R0(J+1)=BLAX(J)
          ELSE IF(II.EQ.LUP+1) THEN
coupling to H-groundstate
            R1(J+1)=BLAX(J)
          ELSE IF(II.EQ.LUP+2) THEN
coupling to Q_EXT
            R_EXT(J+1)=BLAX(J)
          END IF
 3010   CONTINUE
 3001 CONTINUE


      RETURN
      END

C***********************************************************************
      SUBROUTINE IONREC(C,S,SAHA,A,ALPHA,BETA,R0,R1,DENSEL,LUP,LIM,
     &                  F,R_EXT,Q_EXT,
     &                  ALPCR,SCR,SCR_EXT)
C
C     EFFECTIVE IONIZATION AND RECOMBINATION RATE COEFFICIENTS
C     FOR ATOMIC HYDROGEN
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(40,40),S(40),SAHA(40),A(40,40),ALPHA(40),BETA(40),
     &          R0(40),R1(40),R_EXT(40),Q_EXT(40),F(40,40)
C
      DO 5000 I=LUP+1,LIM
        R0(I)=1.0*SAHA(I)
        R1(I)=0.0
        R_EXT(I)=0.0
 5000 CONTINUE
C
      SCR=S(1)
      DO 5001 I=2,LUP
        SUSCR=C(1,I)-R1(I)*(F(I,1)*DENSEL+A(I,1))
 5001 SCR=SCR+SUSCR

      SCR_EXT=0.00
      DO 5002 I=2,LUP
        SUSRAD=Q_EXT(I)-R_EXT(I)*(F(I,1)*DENSEL+A(I,1))
 5002 SCR_EXT=SCR_EXT+SUSRAD

      ALPCR1=DENSEL*ALPHA(1)+BETA(1)

      ALPCR2=0.0
      DO 5003 I=2,LIM

 5003 ALPCR2=ALPCR2+R0(I)*(DENSEL*F(I,1)+A(I,1))

      ALPCR=ALPCR1+ALPCR2
      RETURN
      END

C***********************************************************************
      SUBROUTINE E_IONREC(C,S,SAHA,A,ALPHA,BETA,R0,R1,DENSEL,LUP,LIM,
     &                    F,R_EXT,Q_EXT,E_AT,
     &                    ALPCR,      SCR,    SCR_EXT,
     &                    E_ALPCR,  E_SCR,  E_SCR_EXT,
     &                    E_ALPCR_T,E_SCR_T,E_SCR_EXT_T)
C
C     EFFECTIVE ELECTRON ENERGY LOSS IONIZATION AND RECOMBINATION
C     RATE COEFFICIENTS FOR ATOMIC HYDROGEN
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(40,40),S(40),SAHA(40),A(40,40),ALPHA(40),BETA(40),
     &          R0(40),R1(40),R_EXT(40),Q_EXT(40),F(40,40),E_AT(40)
C
      DO 5000 I=LUP+1,LIM
        R0(I)=1.0*SAHA(I)
        R1(I)=0.0
        R_EXT(I)=0.0
 5000 CONTINUE
C
      E_SCR=0.0
      E_SCR_EXT=0.0
      E_ALPCR=0.0
C
C  FOR TESTING
      E_SCR_T=0.0
      E_SCR_EXT_T=0.0
      E_ALPCR_T=0.0
C
C  EFFECTIV ELECTRON COOLING CORRESPONDING TO
C  "ORDINARY" COUPLING TO GROUND STATE S(I)
C
C  PART I
C  1--> inf.  R1(1)=1
      UH=13.595
      DE=(E_AT(1)-UH)
      E_SCR  =S(1)*DE
C  PART II
C  i1--> inf.  R1(i1)=...
      do i1=2,lup
        DE=(E_AT(i1)-UH)
        E_SCR  =E_SCR+r1(i1)*densel*S(I1)*DE
      enddo
C  PART III
C  1--> I2.  R1(1)=1
      DO I2=2,LUP
        DE=(E_AT(1)-E_AT(I2))
        SUSCR=C(1,I2)*DE
        E_SCR=E_SCR+SUSCR
      ENDDO
C  PART IV
C  i1--> i2. i1<i2,  R1(i1)=...
      DO  I1=2,LUP-1
        DO  I2=I1+1,LUP
          DE=(E_AT(I1)-E_AT(I2))
          SUSCR=r1(i1)*densel*C(I1,I2)*DE
          E_SCR=E_SCR+SUSCR
        ENDDO
      ENDDO
C  PART V
C  i2--> i1.  i2>i1, R1(i2)=...
C            (includes i1=1, i.e., invers to PART III)
      DO I1=1,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I2)-E_AT(I1))
          SUSCR  =r1(i2)*densel*F(I2,I1)*DE
          E_SCR=E_SCR+SUSCR
C  separte treatment of radiation losses alone, only for testing
          SUSCR_T=r1(i2)*A(I2,I1)*(-1.)*DE
          E_SCR_T=E_SCR_T+SUSCR_T
        ENDDO
      ENDDO

C  for test only.  evaluate second formular for E_SCR,
C                  using radiation loss E_SCR_T and effective rate SCR
      UH=13.595
      DE=(E_AT(1)-UH)
      E_SCR_T=E_SCR_T+SCR *DE

C  contribution for "ordinary" coupling to ground state done

c  next: contribution for coupling to external source Q_EXT
c        use same codes as for E_SCR, but with Q_EXT(K) <-- C(1,K)
c                                     and S(1)=0
c        and, of course, R_EXT(k) instead of R1(k)
C        furthermore: no electron energy losses associated with Q_EXT: E_Q_EXT=0
C
C
C  PART I
C  zero, because S(1)=0
C  PART II
C  i1--> inf.  R_EXT(i1)=...
      do i1=2,lup
        DE=(E_AT(i1)-UH)
        E_SCR_EXT  =E_SCR_EXT+R_EXT(i1)*densel*S(I1)*DE
      enddo
C  PART III
C  1--> I2.  R_EXT(1)=1
      DO I2=2,LUP
C  NO ELECTRON ENERGY LOSS IS ASSOCIATED WITH Q_EXT
C  STRICTLY: TOGETHER WITH Q_EXT THERE SHOULD ALSO BE AN E_Q_EXT
C  FOR EXTERNAL SOURCES.
        DE=0
        SUSCR=Q_EXT(I2)*DE
        E_SCR_EXT=E_SCR_EXT+SUSCR
      ENDDO
C  PART IV
C  i1--> i2. i1<i2,  R_EXT(i1)=...
      DO I1=2,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I1)-E_AT(I2))
          SUSCR=R_EXT(i1)*densel*C(I1,I2)*DE
          E_SCR_EXT=E_SCR_EXT+SUSCR
        ENDDO
      ENDDO
C  PART V
C  i2--> i1.  i2>i1, R_EXT(i2)=...
C            (includes i1=1, i.e., invers to PART III)
      DO I1=1,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I2)-E_AT(I1))
          SUSCR  =R_EXT(i2)*densel*F(I2,I1)*DE
          E_SCR_EXT=E_SCR_EXT+SUSCR
C  separte treatment of radiation losses alone, only for testing
C         SUSCR_T=R_EXT(i2)*A(I2,I1)*(-1.)*DE
C         E_SCR_EXT_T=E_SCR_EXT_T+SUSCR_T
        ENDDO
      ENDDO

C  for test only.  evaluate second formular for E_SCR_EXT,
C                  using radiation loss E_SCR_EXT_T and effective rate SCR_EXT
C     UH=13.595
C     DE=(E_AT(1)-UH)
C     E_SCR_EXT_T=E_SCR_EXT_T+SCR_EXT *DE

c  this test only works if in PART III one would set:
C       DE=(E_AT(1)-E_AT(I2))
C  But for external sources this is not necessarily the case
C  test done, 11.01.05,  o.k.,  then test switched off.

c  next: contribution for coupling to H+

C  TO BE WRITTEN  E_ALPRC=....

      RETURN
      END

C********************************************************************
C            C(I,J)  FROM  JOHNSON
C *******************************************************************
      SUBROUTINE JOHN(TEMP,CJ,OSC)
C
C
C     C(I,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN I -> J
C
C     L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972).
C
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 OSC(40,40),CJ(40,40)
      TE=TEMP*1.1605E4
C
      DO 1 I=1,40
      DO 2 J=1,40
      P=I
      Q=J
      IF (I.EQ.1) THEN
      BN=-0.603
      ELSE
      BN=(4.0-18.63/P+36.24/P**2-28.09/P**3)/P
      END IF
      EP=13.6/P**2
      EQ=13.6/Q**2
      U=(EP-EQ)/TEMP
C
      IF (U.GT.0) THEN
      X=1-P**2/Q**2
      Y=U
      R=1.94/P**1.57*X
      Z=R+Y
      A=2*P**2*OSC(I,J)/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X+BN/X**2)
      C1=2*P**2/X
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z
      CALL EXPI(Y1,E1Y,ICON)
      CALL EXPI(Z1,E1Z,ICON)
      E1Y=-E1Y

      E1Z=-E1Z

      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=(1/Y+0.5)*E1Y-(1/Z+0.5)*E1Z
      E2=E2Y/Y-E2Z/Z

      CJ(I,J)=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)
      ELSE
      CJ(I,J)=0.0

      END IF
C
    2 CONTINUE
    1 CONTINUE
C
C
      RETURN
      END
c
      subroutine LAX(A,N1,N,B,eps,ifl,is,vw,ip,icon)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      double precision a(n1,n1),B(*),vw(*),d
      dimension ip(*)
      dimension iw(100),indx(100)
      if (n1.gt.100) then
        write (iunout,*) 'error in lax'
      stop
      endif
      call galpd(a,n1,n,b,iw,ier)
      if (ier.eq.1) then
        write (iunout,*) 'error in lax, matrix ist singulaer'
      endif
      return
      end
c
      SUBROUTINE GALPD(A,NA,NG,B,IW,IER)
C
C***********************************************************************
C*   GAUSS-ALGORITHMUS ZUR LOESUNG LINEARER GLEICHUNGS-SYSTEME MIT     *
C*   PIVOTIERUNG.                                                      *
C*   GENAUIGKEIT:   DOUBLD-PRECISION                   (01.07.1991)    *
C***********************************************************************
C    A(NA,NA): KOEFFIZIENTEN-MATRIX DES GLEICHUNGS-SYSTEMS
C    NA      : DIMENSION VON A WIE IM AUFRUFENDEN PROGRAMM ANGEGEBEN
C    NG      : ANZAHL DER UNBEKANNTEN   (NG <= NA)
C    B(NG)   : ELEMENTE DER RECHTEN SEITE DES GLEICHUNGS-SYSTEMS
C    IW(NG)  : INTEGER-HILFS-ARRAY FUER EINE MOEGLICHE PROGRAMM-
C              INTERNE UMNUMERIERUNG DER GLEICHUNGEN
C    IER     : ERROR-INDEX (IER = 1: MATRIX SINGULAER)
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NA,NA),B(NG),IW(NG)
      DATA ZERO /1.D-71/
      IER=0
C
C     ******************************************************************
C     DER FALL:   NG = 2
C     ******************************************************************
C
      IF(NG.EQ.2) THEN
                  H=A(1,1)*A(2,2)-A(2,1)*A(1,2)
                  AI=B(1)*A(2,2)-B(2)*A(1,2)
                  AK=A(1,1)*B(2)-A(2,1)*B(1)
                  B(1)=AI/H
                  B(2)=AK/H
                  RETURN
                  ENDIF
C
C     ******************************************************************
C     NUMMERN DER UNBEKANNTEN AUF IW ABSPEICHERN.
C     ******************************************************************
C
      DO 1 K=1,NG
    1    IW(K)=K
C
C     ******************************************************************
C     DIE A-MATRIX AUF DREIECKS-FORM BRINGEN.
C     ******************************************************************
C
      DO 10 I=1,NG
C
C        ===============================================================
C        Pivot-Element suchen  (  Zeile IZ,  Spalte KS )
C        ===============================================================
C
         AP=0
         IZ=0
         KS=0
         DO 3 M=I,NG
            DO 2 N=I,NG
               AMN=DABS(A(M,N))
               IF(AMN.GT.AP) THEN
                             AP=AMN

                             IZ=M

                             KS=N

                             ENDIF
    2          CONTINUE
    3      CONTINUE
C
C        ============================
C        Zeilen umordnen, wenn IZ > I
C        ============================
C
         IF(IZ.GT.I) THEN
                     DO 4 N=I,NG
                        H=A(I,N)
                        A(I,N)=A(IZ,N)
    4                   A(IZ,N)=H
                     H=B(I)
                     B(I)=B(IZ)
                     B(IZ)=H
                     ENDIF
C
C        ===============================================
C        Spalten umordnen und Unbekannte neu numerieren
C        ===============================================
C
         IF(KS.GT.I) THEN
                     IH=IW(I)
                     IW(I)=IW(KS)
                     IW(KS)=IH
                     DO 5 M=1,NG
                        AK=A(M,KS)
                        A(M,KS)=A(M,I)
    5                   A(M,I)=AK
                     ENDIF
C
C        ===============================================================
C        Total-Pivotierung. Spalte i,  Zeile 1 ... i-1, i+1 ... NG
C        zu Null machen.
C        ===============================================================
C
         AP=A(I,I)
         IF(DABS(AP).LT.ZERO) GOTO 15
         AP=1/AP
         DO 8 M=1,NG
            IF(M.EQ.I) GOTO 8
            IF(DABS(A(M,I)).GT.ZERO) THEN
                                     Q=A(M,I)*AP
                                     DO 7 N=I,NG
    7                                   A(M,N)=A(M,N)-A(I,N)*Q
                                     B(M)=B(M)-B(I)*Q
                                     ENDIF
    8       CONTINUE
   10    CONTINUE
C
C     ******************************************************************
C     WENN A(N,N) ^= 0, KOENNEN DIE UNBEKANNTEN BERECHNET WERDEN.
C     -:  SIE WERDEN ZUNAECHST AUF A(N,I), I=1,...,N GESETZT.
C     -:  DANN DIE BERECHNETEN UNBEKANNTEN IN DER RICHTIGEN ANORDNUNG
C         AUF B(I) SCHREIBEN UND AN DAS AUFRUFENDE PROGRAMM ZURUECKGEBEN
C     ******************************************************************
C
      DO 12 M=1,NG
   12    A(M,M)=B(M)/A(M,M)
      DO 14 M=1,NG
         II=IW(M)
   14    B(II)=A(M,M)
C
       RETURN

C     ******************************************************************
C     MATRIX IST SINGULAER.
C       -: IER = 1 SETZTEN
C       -: RETURN
C     ******************************************************************
C
   15 IER=1
      RETURN
      END


      subroutine expi(x,ei,icon)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c  exponential integral
c  -int exp(-x)/x, von -x nach unendlich     x<0,  identisch mit
c  +int exp(x)/x,  von -unendl. bis x
c
c   PV +int exp(-x)/x von unendl bis -x     x>0 ,  identisch mit
c   PV +int exp(x)/x von -unendl bis x
      double precision mmdei,dei,xx
cdr   if (x.gt.0) then
        xx=x
        dei=-mmdei(xx,ier)
        icon=ier
        ei=dei
cdr   else
cdr     write (*,*) 'argument in expi lt.0, call exit'
cdr     stop
cdr   endif
      return
      end
c
      subroutine aqc8(a,b,f,epsa,epsr,nmin,nmax,S,err,n,icon)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c  S=integral von a bis b, der function f(x) (external).
c  epsa,epsr : absolute and relative errors, input
c   nmin,nmax  min u max anzahl der functionsaufrufe
c   (nmax<511, nmin>15)
c
c output
c S
c err: estim absolut error
c n  anzahl der functionsaufrufe
c icon: error code
      external f
      external midpnt
      call qromo(f,a,b,s,midpnt,epsa)
      return
      end



c
      DOUBLE PRECISION FUNCTION MMDEI (S,IER)
C  exponential integral function
c  imsl routine, dort: mmdei(ipot=2,s,ier), also:
c  s muss gt.0, mmdei ist dann: integral (s bis unendlich) von
c               exp(-t)/t dt
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      X=S
      Y=ABS(X)
      Z=0.25D+0*Y
      IF(Z-1.0D0)11,11,12
   11 VALUE=((((((((((((((((((((-.483702D-8*Z+.2685377D-7)*Z-.11703642D-
     1 6)*Z+.585911692D-6)*Z-.2843937873D-5)*Z+.1284394756D-4)*Z-.547380
     1 648948D-4)*Z+.21992775413732D-3)*Z-.8290046678016D-3)*Z+.00291878
     1 85699858)*Z-.0095523774824170)*Z+.028895941356815)*Z-.08026651003
     1 2735)*Z+.20317460364863)*Z-.46439909294939)*Z+.948148
     1 1480886)*Z-1.706666666669
     1          )*Z+2.6666666666702)*Z-3.5555555555546)*Z+3.999999999999
     1 4)*Z-3.9999999999996)*Z+.57721566490143+LOG(Y)
      VALUE=-VALUE
      GO TO 13
   12 Z=1.0D0/Z
      VALUE=(((((((((((((((((((-.77769007188383D-3*Z+.84295295893998D-2
     1 )*Z-.04272827083185)*Z+.13473041266261)*Z-.29671551148)*Z+.486188
     1 06480858)*Z-.61743468824936)*Z+.62656052544291)*Z-.52203502518727
     1 )*Z+.36771959879483)*Z-.22727998629908)*Z+.12965447884319)*Z
     1 -.72886262704894D-1)*Z+.043418637381012)*Z-.29244589614262D-1
     1 )*Z+.23433749581188D-1)*Z-.023437317333964)*Z+.03124999448124
     1 )*Z-.062499999910765)*Z+.24999999999935)*Z-.20933859981036D-14
c     if (y.gt.160.D0) then
c       value=0
c     else
        VALUE=VALUE*EXP(-Y)
c     endif
   13 VAL=VALUE
      MMDEI=VAL
      ier=0
      RETURN
      END

      SUBROUTINE QROMO(FUNC,A,B,SS,CHOOSE,eps)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (JMAX=14,JMAXP=JMAX+1,KM=4,K=KM+1)
      DIMENSION S(JMAXP),H(JMAXP)
      external choose,func
      H(1)=1.0D0
      DO 11 J=1,JMAX

        CALL CHOOSE(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.0d0,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) then
            RETURN
          endif
        ENDIF
        S(J+1)=S(J)
        H(J+1)=H(J)/9.
11    CONTINUE
      PAUSE 'Too many steps.'
      write(iunout,*) '(W) Too many steps.'
      END

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N

        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I

          DIF=DIFT

        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
ctk       IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)


        ELSE
          DY=D(NS)


          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
C
      SUBROUTINE MIDPNT(FUNC,A,B,S,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      external func
      save
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5*(A+B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
      RETURN
      END
c
