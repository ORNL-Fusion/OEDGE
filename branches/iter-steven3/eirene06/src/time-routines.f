C EIRENE06 COMPILATION
C ===== SOURCE: timea.f
C
C
      SUBROUTINE TIMEA
C
C   1 ST INTERSECTION OF THE RAY X+T*VX,Y+T*VY,Z+T*VZ WITH ONE OF THE
C   ADDITIONAL SURFACES, DEFINED BY 2.ND ORDER EQUATIONS
C   IT IS ALSO CHECKED, WHETHER THIS INTERSECTION TAKES PLACE INSIDE THE
C   SPECIFIED BOUNDARIES OF THOSE SURFACES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CLGIN
      USE CTRIG
      USE COMSPL
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
C
      REAL(DP) :: XB(3), XC(3), XD(3), XE(3), XF(3), XG(3)
      REAL(DP) :: TADD, TMIN, XXR, V, VSAVE, A1, A2, PPP, ROT, VX, VY,
     .          VZ, VXS, VYS, VZS,VVS, TS, T, X, Y, Z, A3, VXJ, VYJ,
     .          VZJ, XLS2, XLS3, XMS3, XMS2, XX, YY, ZZ, XR, YR, ZR,
     .          XN, YN, ZN, TMA, TMI, WR, TST, XMS1, XLS1, TMX, TUP,
     .          DZ3, TL, XS, DZ, TEST4, TEST5, DY, DL, DX1, DX, SG,
     .          VV, VXX, VYY, VZZ, TMT, B, XK5, V445, XN5, V335, XK4,
     .          V334, V3V45, XN4, YS, ZS, AT, XNORM, V224, CG1, CG2,
     .          CG3, AR1, C, D, E, F, G, XK3, HELP, V2V34, XN3, V223
     .          AR3, AR2, V113, V1V23, V223, AR3, DY2
      INTEGER :: NN, NLLLI, NTNEW, NNTCLS, NNR, ICOUNT, IPLGN, ITRII,
     .           ISTS, I1000, J, I, IPERID, NTCELL, NCELL, MSURF,
     .           NLI, NLE, MASURF, ISPZ, NNTCL, LM2, LM1, IAB, JUM
      INTEGER, EXTERNAL :: IDEZ
      LOGICAL :: LGJ, NLTRC, BITGET, LCNDEXP, LTSTCXP
c slmod begin - gfortran
      LOGICAL, ALLOCATABLE :: LMTSRF(:)
c
c     LOGICAL :: LMTSRF(NLIMPS)
c slmod end
      SAVE
C
      ENTRY TIMEA0
c slmod begin - gfortran
      IF (.NOT.ALLOCATED(LMTSRF)) ALLOCATE(LMTSRF(NLIMPS))
c slmod end
C     
      IF (NLIMI.LT.1) RETURN
C
C
      DO 2 J=1,NLIMI
        IF (IGJUM0(J).NE.0) THEN
          IF (NLIMPB >= NLIMPS) THEN
            DO 1 I=0,NLIMI
              IGJUM1(I,J)=1
1           CONTINUE
          ELSE
            DO I=0,NLIMI
              CALL BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
            END DO
          END IF
        ENDIF
C
        ISWICH(1,J)=IDEZ(ILSWCH(J),1,6)
        IF (ISWICH(1,J).EQ.1) ISWICH(1,J)=-1
        IF (ISWICH(1,J).EQ.2) ISWICH(1,J)=1
        ISWICH(2,J)=IDEZ(ILSWCH(J),2,6)
        IF (ISWICH(2,J).EQ.1) ISWICH(2,J)=-1
        IF (ISWICH(2,J).EQ.2) ISWICH(2,J)=1
        ISWICH(3,J)=IDEZ(ILSWCH(J),3,6)
        IF (ISWICH(3,J).EQ.1) ISWICH(3,J)=-1
        IF (ISWICH(3,J).EQ.2) ISWICH(3,J)=1
        ISWICH(4,J)=IDEZ(ILSWCH(J),4,6)
        IF (ISWICH(4,J).EQ.1) ISWICH(4,J)=-1
        IF (ISWICH(4,J).EQ.2) ISWICH(4,J)=1
        ISWICH(5,J)=IDEZ(ILSWCH(J),5,6)
        IF (ISWICH(5,J).EQ.1) ISWICH(5,J)=-1
        IF (ISWICH(5,J).EQ.2) ISWICH(5,J)=1
        ISWICH(6,J)=IDEZ(ILSWCH(J),6,6)
        IF (ISWICH(6,J).EQ.1) ISWICH(6,J)=-1
        IF (ISWICH(6,J).EQ.2) ISWICH(6,J)=1
C
        IF (ISWICH(4,J).NE.0.OR.ISWICH(5,J).NE.0.OR.ISWICH(6,J).NE.0)
     .  THEN
          ILBLCK(J)=IDEZ(ILCELL(J),4,4)
          I1000=1000*ILBLCK(J)
          ILACLL(J)=ILCELL(J)-I1000
        ENDIF
2     CONTINUE
C
      DO 97 J=1,NLIMI
        IF (IGJUM0(J).NE.0) THEN
          IF (LEVGEO.EQ.4) THEN
            DO ITRII=1,NTRII
            DO IPLGN=1,3
              ISTS=ABS(INMTI(IPLGN,ITRII))
              IF (J.EQ.ISTS) GOTO 85
            ENDDO
            ENDDO
          ENDIF
          GOTO 97
        ENDIF
85      IF (RLB(J).LT.2.0) THEN
C
C   SURFACE COEFFICIENTS ARE INPUT
C
C  INDICATE INDEPENDENCE OF SURFACE EQUATION FROM X, Y, Z, RESP.
          IF (A1LM(J).EQ.0..AND.A4LM(J).EQ.0..AND.
     .        A7LM(J).EQ.0..AND.A8LM(J).EQ.0.)
     .    P3(1,J)=1.D55
          IF (A2LM(J).EQ.0..AND.A5LM(J).EQ.0..AND.
     .        A7LM(J).EQ.0..AND.A9LM(J).EQ.0.)
     .    P3(2,J)=1.D55
          IF (A3LM(J).EQ.0..AND.A6LM(J).EQ.0..AND.
     .        A8LM(J).EQ.0..AND.A9LM(J).EQ.0.)
     .    P3(3,J)=1.D55
          GOTO 90
        ELSEIF (RLB(J).GE.2) THEN
C
C  PLANE SURFACE, VERTICES ARE INPUT, SET SURFACE COEFFICIENTS:
C
          A4LM(J)=0.
          A5LM(J)=0.
          A6LM(J)=0.
          A7LM(J)=0.
          A8LM(J)=0.
          A9LM(J)=0.
        ENDIF
C
        IF (RLB(J).GE.2.6.AND.RLB(J).LT.2.9) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE X CO-ORDINATE
          IF (RLB(J).LT.2.75) RLB(J)=1.
          IF (RLB(J).GE.2.75) RLB(J)=1.5
          XLIMS1(1,J)=P1(1,J)
          XLIMS2(1,J)=P2(1,J)
          DZ3=(P1(3,J)-P2(3,J))
          DY2=(P1(2,J)-P2(2,J))
!pb  avoid a0=a1=a2=a3=0
          IF (ABS(DZ3)+ABS(DY2) .LT. EPS12) GOTO 98
          A0LM(J)=DZ3*P1(2,J)-DY2*P1(3,J)
          A1LM(J)=0.
          A2LM(J)=-DZ3
          A3LM(J)=DY2
          YLIMS1(1,J)=MIN(P1(2,J),P2(2,J))
          YLIMS2(1,J)=MAX(P1(2,J),P2(2,J))
          ZLIMS1(1,J)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(1,J)=MAX(P1(3,J),P2(3,J))
!
!          write (iunout,'(a,6es12.4)') 'xylims1,2 ',
!     .      xlims1(1,j),ylims1(1,j),zlims1(1,j),
!     .      xlims2(1,j),ylims2(1,j),zlims2(1,j)
!
          IF (P1(2,J).EQ.P2(2,J)) THEN
            YLIMS1(1,J)=YLIMS1(1,J)-0.1
            YLIMS2(1,J)=YLIMS2(1,J)+0.1
          ENDIF
          IF (P1(3,J).EQ.P2(3,J)) THEN
            ZLIMS1(1,J)=ZLIMS1(1,J)-0.1
            ZLIMS2(1,J)=ZLIMS2(1,J)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND X-INDEPENDENCE OF SURFACE EQUATION
          P3(1,J)=1.D55
          DL=SQRT(DY2**2+DZ3**2)
          DX=XLIMS2(1,J)-XLIMS1(1,J)
          IF (DX.GT.1.D20) DX=XDF
          SAREA(J)=DL*DX
          GOTO 90
        ELSEIF (RLB(J).GE.2.3.AND.RLB(J).LT.2.6) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE Y CO-ORDINATE
          IF (RLB(J).LT.2.45) RLB(J)=1.
          IF (RLB(J).GE.2.45) RLB(J)=1.5
          YLIMS1(1,J)=P1(2,J)
          YLIMS2(1,J)=P2(2,J)
          DZ3=(P1(3,J)-P2(3,J))
          DX1=(P1(1,J)-P2(1,J))
!pb  avoid a0=a1=a2=a3=0
          IF (ABS(DZ3)+ABS(DX1) .LT. EPS12) GOTO 98
          A0LM(J)=DZ3*P1(1,J)-DX1*P1(3,J)
          A1LM(J)=-DZ3
          A2LM(J)=0.
          A3LM(J)=DX1
          XLIMS1(1,J)=MIN(P1(1,J),P2(1,J))
          XLIMS2(1,J)=MAX(P1(1,J),P2(1,J))
          ZLIMS1(1,J)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(1,J)=MAX(P1(3,J),P2(3,J))
!
!          write (iunout,'(a,6es12.4)') 'xylims1,2 ',
!     .      xlims1(1,j),ylims1(1,j),zlims1(1,j),
!     .      xlims2(1,j),ylims2(1,j),zlims2(1,j)
!
          IF (P1(1,J).EQ.P2(1,J)) THEN
            XLIMS1(1,J)=XLIMS1(1,J)-0.1
            XLIMS2(1,J)=XLIMS2(1,J)+0.1
          ENDIF
          IF (P1(3,J).EQ.P2(3,J)) THEN
            ZLIMS1(1,J)=ZLIMS1(1,J)-0.1
            ZLIMS2(1,J)=ZLIMS2(1,J)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND Y-INDEPENDENCE OF SURFACE EQUATION
          P3(2,J)=1.D55
          DL=SQRT(DX1**2+DZ3**2)
          DY=YLIMS2(1,J)-YLIMS1(1,J)
          IF (DY.GT.1.E20) DY=YDF
          SAREA(J)=DL*DY
          GOTO 90
        ELSEIF (RLB(J).GE.2.0.AND.RLB(J).LT.2.3) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE Z CO-ORDINATE
          IF (RLB(J).LT.2.15) RLB(J)=1.
          IF (RLB(J).GE.2.15) RLB(J)=1.5
          ZLIMS1(1,J)=P1(3,J)
          ZLIMS2(1,J)=P2(3,J)
          DY2=(P1(2,J)-P2(2,J))
          DX1=(P1(1,J)-P2(1,J))
!pb  avoid a0=a1=a2=a3=0
          IF (ABS(DY2)+ABS(DX1) .LT. EPS12) GOTO 98
          A0LM(J)=DY2*P1(1,J)-DX1*P1(2,J)
          A1LM(J)=-DY2
          A2LM(J)=DX1
          A3LM(J)=0.
          XLIMS1(1,J)=MIN(P1(1,J),P2(1,J))
          XLIMS2(1,J)=MAX(P1(1,J),P2(1,J))
          YLIMS1(1,J)=MIN(P1(2,J),P2(2,J))
          YLIMS2(1,J)=MAX(P1(2,J),P2(2,J))
!
!          write (iunout,'(a,6es12.4)') 'xylims1,2 ',
!     .      xlims1(1,j),ylims1(1,j),zlims1(1,j),
!     .      xlims2(1,j),ylims2(1,j),zlims2(1,j)
!
          IF (P1(1,J).EQ.P2(1,J)) THEN
!pb          IF (xlims1(1,J).EQ.xlims2(1,J)) THEN
            XLIMS1(1,J)=XLIMS1(1,J)-0.1
            XLIMS2(1,J)=XLIMS2(1,J)+0.1
          ENDIF
          IF (P1(2,J).EQ.P2(2,J)) THEN
!pb          IF (ylims1(1,J).EQ.ylims2(1,J)) THEN
            YLIMS1(1,J)=YLIMS1(1,J)-0.1
            YLIMS2(1,J)=YLIMS2(1,J)+0.1
          ENDIF

C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND Z-INDEPENDENCE OF SURFACE EQUATION
          P3(3,J)=1.D55
          DL=SQRT(DX1**2+DY2**2)
          IF (NLTRZ) THEN
            DZ=ZLIMS2(1,J)-ZLIMS1(1,J)
            IF (DZ.GT.1.D20) DZ=ZDF
          ELSEIF (NLTRA) THEN
            XS=(P1(1,J)+P2(1,J))*0.5+RMTOR
            IF (ILTOR(J).GT.0) THEN
              DZ=ZLIMS2(1,J)-ZLIMS1(1,J)
C             IF (DZ.GT.1.D20) ????
            ELSEIF (ILTOR(J).EQ.0) THEN
              DZ=ZLIMS2(1,J)-ZLIMS1(1,J)
              IF (DZ.GT.1.D20) DZ=XS*TANAL/ALPHA*PI2A
            ENDIF
          ENDIF
          SAREA(J)=DL*DZ
          GOTO 90
        ELSEIF (RLB(J).GE.3.AND.RLB(J).LT.4.) THEN
C  1 TRIANGLE
          P4(1,J)=P3(1,J)
          P4(2,J)=P3(2,J)
          P4(3,J)=P3(3,J)
          P5(1,J)=P4(1,J)
          P5(2,J)=P4(2,J)
          P5(3,J)=P4(3,J)
        ELSEIF (RLB(J).GE.4..AND.RLB(J).LT.5) THEN
C  1 QUADRANGLE = 2 TRIANGLES
          P5(1,J)=P4(1,J)
          P5(2,J)=P4(2,J)
          P5(3,J)=P4(3,J)
        ENDIF
C
        TEST4=0.
        TEST5=0.
C
C  SET PLANE SURFACE FROM CO-ORDINATES OF THE FIRST 3 VERTICES
C  RLB.GE.3.0 AT THIS POINT
C
        DO 82 I=1,3
          XB(I)=P1(I,J)-P3(I,J)
          XC(I)=P2(I,J)-P3(I,J)
          XD(I)=P2(I,J)-P4(I,J)
          XE(I)=P3(I,J)-P4(I,J)
          XF(I)=P3(I,J)-P5(I,J)
          XG(I)=P4(I,J)-P5(I,J)
82      CONTINUE
C
        A1LM(J)=XB(2)*XC(3)-XB(3)*XC(2)
        A2LM(J)=XB(3)*XC(1)-XB(1)*XC(3)
        A3LM(J)=XB(1)*XC(2)-XB(2)*XC(1)
        A0LM(J)=-(A1LM(J)*P1(1,J)+A2LM(J)*P1(2,J)+A3LM(J)*P1(3,J))
C
C   SURFACE AREA: SAREA
        B=SQRT(XB(1)**2+XB(2)**2+XB(3)**2+EPS60)
        C=SQRT(XC(1)**2+XC(2)**2+XC(3)**2+EPS60)
        D=SQRT(XD(1)**2+XD(2)**2+XD(3)**2+EPS60)
        E=SQRT(XE(1)**2+XE(2)**2+XE(3)**2+EPS60)
        F=SQRT(XF(1)**2+XF(2)**2+XF(3)**2+EPS60)
        G=SQRT(XG(1)**2+XG(2)**2+XG(3)**2+EPS60)
        CG1=(XB(1)*XC(1)+XB(2)*XC(2)+XB(3)*XC(3))/B/C
        CG2=(XD(1)*XE(1)+XD(2)*XE(2)+XD(3)*XE(3))/D/E
        CG3=(XF(1)*XG(1)+XF(2)*XG(2)+XF(3)*XG(3))/F/G
        AR1=0.5*B*C*SQRT(1.-CG1*CG1+EPS60)
        AR2=0.5*D*E*SQRT(1.-CG2*CG2+EPS60)
        AR3=0.5*F*G*SQRT(1.-CG3*CG3+EPS60)
        SAREA(J)=AR1+AR2+AR3
C
C   TEST, IF 4TH  AND 5TH POINT ON SURFACE
        TEST4=A0LM(J)+A1LM(J)*P4(1,J)+A2LM(J)*P4(2,J)+A3LM(J)*P4(3,J)
        TEST5=A0LM(J)+A1LM(J)*P5(1,J)+A2LM(J)*P5(2,J)+A3LM(J)*P5(3,J)
C
C  NEW CO-ORDINATE SYSTEM IN P3,P1,P2 ; P4,P2,P3 ; P5,P3,P4
C    XS = P3 + XLS1*(P1-P3) + XMS1*(P2-P3)
C    XS = P4 + XLS2*(P2-P4) + XMS2*(P3-P4)
C    XS = P5 + XLS3*(P3-P5) + XMS3*(P4-P5)
C  PREPARE ARRAYS FOR COMPUTATION OF XLS1,...XMS3
        V1V23=XB(1)*XC(1)+XB(2)*XC(2)+XB(3)*XC(3)
        V113=XB(1)*XB(1)+XB(2)*XB(2)+XB(3)*XB(3)
        V223=XC(1)*XC(1)+XC(2)*XC(2)+XC(3)*XC(3)
C
        IF (ABS(V113).LE.EPS12) GOTO 98
        HELP=V1V23*V1V23/V113-V223
        IF (ABS(HELP).LE.EPS12) GOTO 98
        XK3=V1V23/V113
        XN3=1./HELP
        PS13(1,J)=(XB(1)*XK3-XC(1))*XN3
        PS13(2,J)=(XB(2)*XK3-XC(2))*XN3
        PS13(3,J)=(XB(3)*XK3-XC(3))*XN3
        P1A(J)=-P3(1,J)*PS13(1,J)-P3(2,J)*PS13(2,J)-P3(3,J)*PS13(3,J)
C
        IF (ABS(V223).LE.EPS12) GOTO 98
        HELP=V1V23*V1V23/V223-V113
        IF (ABS(HELP).LE.EPS12) GOTO 98
        XK3=V1V23/V223
        XN3=1./HELP
        PS23(1,J)=(XC(1)*XK3-XB(1))*XN3
        PS23(2,J)=(XC(2)*XK3-XB(2))*XN3
        PS23(3,J)=(XC(3)*XK3-XB(3))*XN3
        P2A(J)=-P3(1,J)*PS23(1,J)-P3(2,J)*PS23(2,J)-P3(3,J)*PS23(3,J)
C
        IF (RLB(J).GE.4) THEN
          V2V34=XD(1)*XE(1)+XD(2)*XE(2)+XD(3)*XE(3)
          V224=XD(1)*XD(1)+XD(2)*XD(2)+XD(3)*XD(3)
          V334=XE(1)*XE(1)+XE(2)*XE(2)+XE(3)*XE(3)
C
          IF (ABS(V224).LE.EPS12) GOTO 98
          HELP=V2V34*V2V34/V224-V334
          IF (ABS(HELP).LE.EPS12) GOTO 98
          XK4=V2V34/V224
          XN4=1./HELP
          PS24(1,J)=(XD(1)*XK4-XE(1))*XN4
          PS24(2,J)=(XD(2)*XK4-XE(2))*XN4
          PS24(3,J)=(XD(3)*XK4-XE(3))*XN4
          P1B(J)=-P4(1,J)*PS24(1,J)-P4(2,J)*PS24(2,J)-P4(3,J)*PS24(3,J)
C
          IF (ABS(V334).LE.EPS12) GOTO 98
          HELP=V2V34*V2V34/V334-V224
          IF (ABS(HELP).LE.EPS12) GOTO 98
          XK4=V2V34/V334
          XN4=1./HELP
          PS34(1,J)=(XE(1)*XK4-XD(1))*XN4
          PS34(2,J)=(XE(2)*XK4-XD(2))*XN4
          PS34(3,J)=(XE(3)*XK4-XD(3))*XN4
          P2B(J)=-P4(1,J)*PS34(1,J)-P4(2,J)*PS34(2,J)-P4(3,J)*PS34(3,J)
C
          IF (RLB(J).GE.5) THEN
            V3V45=XF(1)*XG(1)+XF(2)*XG(2)+XF(3)*XG(3)
            V335=XF(1)*XF(1)+XF(2)*XF(2)+XF(3)*XF(3)
            V445=XG(1)*XG(1)+XG(2)*XG(2)+XG(3)*XG(3)
C
            IF (ABS(V335).LE.EPS12) GOTO 98
            HELP=V3V45*V3V45/V335-V445
            IF (ABS(HELP).LE.EPS12) GOTO 98
            XK5=V3V45/V335
            XN5=1./HELP
            PS35(1,J)=(XF(1)*XK5-XG(1))*XN5
            PS35(2,J)=(XF(2)*XK5-XG(2))*XN5
            PS35(3,J)=(XF(3)*XK5-XG(3))*XN5
            P1C(J)=-P5(1,J)*PS35(1,J)-P5(2,J)*PS35(2,J)-
     -              P5(3,J)*PS35(3,J)
C
            IF (ABS(V445).LE.EPS12) GOTO 98
            HELP=V3V45*V3V45/V445-V335
            IF (ABS(HELP).LE.EPS12) GOTO 98
            XK5=V3V45/V445
            XN5=1./HELP
            PS45(1,J)=(XG(1)*XK5-XF(1))*XN5
            PS45(2,J)=(XG(2)*XK5-XF(2))*XN5
            PS45(3,J)=(XG(3)*XK5-XF(3))*XN5
            P2C(J)=-P5(1,J)*PS45(1,J)-P5(2,J)*PS45(2,J)-
     -              P5(3,J)*PS45(3,J)
          ENDIF
        ENDIF
C
        IF (ABS(TEST4).GT.EPS10) THEN
          WRITE (iunout,*) 'WARNING FROM TIMEA0, TEST FOR 4.TH POINT:'
          WRITE (iunout,*) 'SF. NUMBER, TEST= ',J,TEST4
        ELSEIF (ABS(TEST5).GT.EPS10) THEN
          WRITE (iunout,*) 'WARNING FROM TIMEA0, TEST FOR 5.TH POINT:'
          WRITE (iunout,*) 'SF. NUMBER, TEST= ',J,TEST5
        ENDIF
C
C  SET SOME MORE ASSISTENT SURFACE DATA AND CHECK CONSISTENCY
C  ST. NO. 90 --- 99
C
90      CONTINUE
C
        RLBNOT(J)=RLB(J).EQ.1.5.OR.RLB(J).EQ.2.5.OR.RLB(J).EQ.3.5.OR.
     .            RLB(J).EQ.4.5.OR.RLB(J).EQ.5.5
C
        IF (ILSWCH(J).NE.0.AND.ILIIN(J).EQ.0) THEN
          WRITE (iunout,*) 'EXIT FROM TIMEA0: SURFACE NO J IS OPERATING'
          WRITE (iunout,*) 'A SWITCH BUT IS NOT SEEN BY HISTORY'
          WRITE (iunout,*) 'J= ',J
          CALL EXIT_OWN(1)
        ENDIF
C
        IF (ILSWCH(J).NE.0.AND.ILIIN(J).GT.0) THEN
          DO ISPZ=1,NSPTOT
            IF (TRANSP(ISPZ,1,J).NE.0.D0.OR.TRANSP(ISPZ,2,J).NE.0.D0)
     .      GOTO 95
          ENDDO
          GOTO 96
95        WRITE (iunout,*) 'EXIT FROM TIMEA0: SURFACE NO J IS OPERATING'
          WRITE (iunout,*) 'A SWITCH BUT IS SOMETIMES TRANSPARENT AND '
          WRITE (iunout,*) 
     .      'SOMETIMES REFLECTING (SEMI-TRANSPARENCY OPTION)'
          WRITE (iunout,*) 
     .      'POSSIBLE FIXES: MANUAL, CHAPTER 2, SECTION 6 '
          WRITE (iunout,*) 'J= ',J
          CALL EXIT_OWN(1)
        ENDIF
C
96      GOTO 94
C
98      WRITE (iunout,*) 
     .    'WARNING FROM TIMEA0, SURFACE NO J IS TURNED OFF'
        WRITE (iunout,*) 'BECAUSE ILL DEFINED '
        WRITE (iunout,*) 'J= ',J
        IGJUM0(J)=1
        IF (NLIMPB >= NLIMPS) THEN
          DO IAB=0,NLIMI
            IGJUM1(IAB,J)=1
          END DO
        ELSE
          DO IAB=0,NLIMI
            CALL BITSET (IGJUM1,0,NLIMPS,IAB,J,1,NBITS)
          END DO
        END IF
        GOTO 97
C
94      CONTINUE
C
        JUMLIM(J)=0
        ALM(J)=2.*A4LM(J)
        BLM(J)=2.*A5LM(J)
        CLM(J)=2.*A6LM(J)
        IF (A4LM(J).EQ.0..AND.A5LM(J).EQ.0..AND.A6LM(J).EQ.0..AND.
     .      A7LM(J).EQ.0..AND.A8LM(J).EQ.0..AND.A9LM(J).EQ.0.) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM1(J,J)=1
          ELSE
            CALL BITSET (IGJUM1,0,NLIMPS,J,J,1,NBITS)
          END IF
          AT=MAX(ABS(A1LM(J)),ABS(A2LM(J)),ABS(A3LM(J)))
          IF (ABS(A1LM(J)).EQ.AT) JUMLIM(J)=1
          IF (ABS(A2LM(J)).EQ.AT) JUMLIM(J)=2
          IF (ABS(A3LM(J)).EQ.AT) JUMLIM(J)=3
          XNORM=SQRT(A1LM(J)*A1LM(J)+A2LM(J)*A2LM(J)+A3LM(J)*A3LM(J))
          A0LM(J)=A0LM(J)/XNORM
          A1LM(J)=A1LM(J)/XNORM
          A2LM(J)=A2LM(J)/XNORM
          A3LM(J)=A3LM(J)/XNORM
          JUM=JUMLIM(J)
          GOTO (91,92,93),JUM
91          ALM(J)=-A0LM(J)/A1LM(J)
            BLM(J)=-A2LM(J)/A1LM(J)
            CLM(J)=-A3LM(J)/A1LM(J)
          GOTO 97
92          ALM(J)=-A0LM(J)/A2LM(J)
            BLM(J)=-A1LM(J)/A2LM(J)
            CLM(J)=-A3LM(J)/A2LM(J)
          GOTO 97
93          ALM(J)=-A0LM(J)/A3LM(J)
            BLM(J)=-A1LM(J)/A3LM(J)
            CLM(J)=-A2LM(J)/A3LM(J)
        ENDIF
97    CONTINUE
C
      CALL LEER(2)
99    CONTINUE
      RETURN
C
      ENTRY TIMEA1(MSURF,NCELL,NLI,NLE,NTCELL,IPERID,XX,YY,ZZ,TMT,
     .             VXX,VYY,VZZ,VV,
     .             MASURF,XR,YR,ZR,SG,TL,NLTRC,LCNDEXP)
c slmod begin - gfortran
      IF (.NOT.ALLOCATED(LMTSRF)) ALLOCATE(LMTSRF(NLIMPS))
c slmod end
      LM1=NLI
      LM2=NLE
      LMTSRF=.FALSE.
      LCNDEXP=.FALSE.
C  PARTICLE ON STANDARD SURFACE?
      IF (MSURF.GT.NLIM) MSURF=0
C  FIND LOCAL COORDINATE SYSTEM IN CASE OF NLTRA
      NNTCL=1
      IF (NLTRA.AND.NLTOR) THEN
        NNTCL=NTCELL
      ELSEIF (NLTRA.AND..NOT.NLTOR) THEN
        NNTCL=IPERID
      ENDIF
C  SAVE INITIAL CO-ORDINATES
      XS=XX
      YS=YY
      ZS=ZZ
      TS=TMT
      VXS=VXX
      VYS=VYY
      VZS=VZZ
      VVS=VV
      NNTCLS=NNTCL
C  WORKING CO-ORDINATES
      X=XX
      Y=YY
      Z=ZZ
      T=TMT
      VX=VXX
      VY=VYY
      VZ=VZZ
      V=VV
      NN=NNTCL
C
      NLLLI=MSURF
      IF (NLTRC) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'TIMEA ,X,Y,Z,T ',X,Y,Z,T
        WRITE (iunout,*) '      VX,VY,VZ,V ',VX,VY,VZ,V
        IF (NLTRA) WRITE (iunout,*) 'MSURF,NNTCL ',MSURF,NNTCL
        IF (.NOT.NLTRA) WRITE (iunout,*) 'MSURF ',MSURF
      ENDIF
C
      TADD=0.
C
1000  TMIN=1.D30
      TL=1.D30
      MASURF=0
C
C  LOOP OVER SURFACE NUMBER, DO 100
C
      DO 100 J=LM1,LM2
        IF (IGJUM0(J).NE.0) GOTO 100
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM1(NLLLI,J) .NE. 0) GOTO 100
        ELSE
          IF (BITGET(IGJUM1,0,NLIMPS,NLLLI,J,NBITS)) GOTO 100
        END IF
        IF (NCELL.LE.NOPTIM) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IF (IGJUM3(NCELL,J).NE.0) GOTO 100
          ELSE
            IF (BITGET(IGJUM3,0,NOPTIM,NCELL,J,NBITS)) GOTO 100
          END IF
        ENDIF
C
C  FOR NLTRA OPTION ONLY:
C  X,Z,VX AND VZ ARE GIVEN IN TOROIDAL CELL NN,
C  TRANSFORM CO-ORDINATES FOR THIS TRACK FROM LOCAL SYSTEM NN
C  TO THE LOCAL SYSTEM ILTOR(J), IN WHICH SURFACE J IS GIVEN
C  IF (ILTOR(J).LE.0) THIS SURFACE HAS TOROIDAL SYMMETRY
C
        IF (NLTRA) THEN
          IF (ILTOR(J).GT.0.AND.ILTOR(J).NE.NN) THEN
            CALL FZRTOR (X,Z,NN,XXR,PPP,NTNEW,.FALSE.,0)
            CALL FZRTRI (X,Z,ILTOR(J),XXR,PPP,NTNEW)
            ROT=2.*(NN-ILTOR(J))*ALPHA
            VSAVE=VX
            VX=COS(ROT)*VSAVE-SIN(ROT)*VZ
            VZ=SIN(ROT)*VSAVE+COS(ROT)*VZ
            NN=ILTOR(J)
C           WRITE (iunout,*) 'J,ILTOR(J),X,Z,VX,VZ ',
C    .                        J,ILTOR(J),X,Z,VX,VZ
C  X,Z,VX AND VZ ARE NOW GIVEN IN CELL NN=ILTOR(J). SO ARE THE COEFFICIENTS
C  OF SURFACE NO. J. FIND INTERSECTION IN THIS LOCAL SYSTEM
          ELSEIF (ILTOR(J).EQ.0) THEN
C  TOROIDALLY SYMMETRIC SURFACE, SURFACE COEFFICIENTS ARE THE SAME
C  IN EACH TOROIDAL CELL, THUS ESPECIALLY IN CELL NN
            NN=NNTCLS
            X=XS
            Z=ZS
            VX=VXS
            VZ=VZS
C           WRITE (iunout,*) 'J,ILTOR(J),X,Z,VX,VZ ',
C    .                        J,ILTOR(J),X,Z,VX,VZ
          ENDIF
        ENDIF
C
C  FIND INTERSECTION TIME TMX WITH BOUNDARY NO. J
C
        LTSTCXP=.FALSE.
        A1=0.
        GOTO (60,63,66),JUMLIM(J)
C  A1*TMX*TMX+A2*TMX+A3=0
        A1=(A4LM(J)*VX+A7LM(J)*VY+A8LM(J)*VZ)*VX+
     .     (A5LM(J)*VY+A9LM(J)*VZ)*VY+A6LM(J)*VZ*VZ
        A2=(A1LM(J)+ALM(J)*X)*VX+(A2LM(J)+BLM(J)*Y)*VY+
     .     (A3LM(J)+CLM(J)*Z)*VZ+
     .      A7LM(J)*(VX*Y+VY*X)+A8LM(J)*(VX*Z+VZ*X)+A9LM(J)*(VY*Z+VZ*Y)
C  A3 =  0. ?
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM2(NLLLI,J).NE.0) GOTO 40
        ELSE
          IF (BITGET(IGJUM2,0,NLIMPS,NLLLI,J,NBITS)) GOTO 40
        END IF
C  NO
        A3=A0LM(J)+(A1LM(J)+A4LM(J)*X+A7LM(J)*Y+A8LM(J)*Z)*X
     .            +(A2LM(J)+A5LM(J)*Y+A9LM(J)*Z)*Y
     .            +(A3LM(J)+A6LM(J)*Z)*Z
        IF (A1.EQ.0.) GOTO 50
        F=-A2/(A1+A1)
        G=F*F-A3/A1
        IF (G.LT.0.) GOTO 100
        G=SQRT(G)
        TMA=F+G
        TMI=F-G
        IF (NLTRC) WRITE (iunout,*) 'TIMEA, J,TMI,TMA ',J,TMI,TMA
!pb        IF (TMA.LE.EPS12.OR.TMI.GT.TMIN) GOTO 100
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMI.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
!  NO SWITCHING SURFACE FOR CONDITIONAL EXPECTATION ESTIMATOR
            GOTO 100
          ELSE
!  SWITCHING SURFACE FOR CONDITIONAL EXPECTATION ESTIMATOR
            LTSTCXP=.TRUE.
          END IF
        END IF
        IF (RLB(J).LT.0.) GOTO 21
        IF (RLB(J).GT.0.) GOTO 31
C
C  RLB(J) .EQ. 0.  STATEMENT 11---20
C
11      IF (TMI.LE.EPS12) GOTO 12
        IF (LTSTCXP) THEN
          LCNDEXP =.TRUE.
          GOTO 100
        END IF
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        TMX=TMI
        GOTO 70
C
12      IF (TMA.GT.TMIN) THEN
          IF (LTSTCXP) LCNDEXP =.TRUE.
          GOTO 100
        END IF
        WR=A2+A1*(TMA+TMA)
16      XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
18      TMX=TMA
        GOTO 70
C
C  CHECK BOUNDARY INEQUALITIES OF SURFACE
C  LGJ=.TRUE.: INTERSECTION POINT IS STILL VALID
C  LGJ=.FALSE.: INTERSECTION POINT IS OUTSIDE THE SPECIFIED AREA
C
C  RLB(J) .LT. 0.  STATEMENT 21---30
C
21      IF (TMI.LE.EPS12) GOTO 22
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        TUP=TMI
        ICOUNT=1
        GOTO 28
C
22      IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2+A1*(TMA+TMA)
26      XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
        TUP=TMA
        ICOUNT=2
C
28      LGJ=.TRUE.
        IF (ILIN(J).GT.0) THEN
          I=0
27        I=1+I
          TST=ALIMS(I,J)+XLIMS(I,J)*XN+YLIMS(I,J)*YN+ZLIMS(I,J)*ZN
          LGJ=TST.LE.0.
          IF (LGJ.AND.I.LT.ILIN(J)) GOTO 27
        ENDIF
        IF (LGJ.AND.ISCN(J).GT.0) THEN
          I=0
29        I=1+I
          TST=ALIMS0(I,J)+
     .        XN*(XLIMS1(I,J)+XN*XLIMS2(I,J)+YN*XLIMS3(I,J))+
     .        YN*(YLIMS1(I,J)+YN*YLIMS2(I,J)+ZN*ZLIMS3(I,J))+
     .        ZN*(ZLIMS1(I,J)+ZN*ZLIMS2(I,J)+XN*YLIMS3(I,J))
          LGJ=TST.LE.0.
          IF (LGJ.AND.I.LT.ISCN(J)) GOTO 29
        ENDIF
        IF (NLPRCS(J).AND.LGJ) LCNDEXP=.TRUE.
        IF (LGJ) THEN
          IF (LTSTCXP) GOTO 100
          TMX=TUP
          GOTO 70
        ELSEIF (ICOUNT.EQ.1) THEN
          GOTO 22
        ENDIF
        GOTO 100
C
C   RLB(J) .GT. 0.  STATEMENT 31---40
C
31      IF (TMI.LE.EPS12) GOTO 32
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        LGJ=XN.LE.XLIMS2(1,J).AND.XN.GE.XLIMS1(1,J).AND.
     .      YN.LE.YLIMS2(1,J).AND.YN.GE.YLIMS1(1,J).AND.
     .      ZN.LE.ZLIMS2(1,J).AND.ZN.GE.ZLIMS1(1,J)
        IF (RLBNOT(J)) LGJ=.NOT.LGJ
        IF (NLPRCS(J).AND.LGJ) LCNDEXP=.TRUE.
        IF (LTSTCXP) GOTO 100
        TMX=TMI
        IF (LGJ) GOTO 70
C
32      IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2+A1*(TMA+TMA)
36      XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
38      CONTINUE
        IF (RLB(J).LT.2.) THEN
          LGJ=XN.LE.XLIMS2(1,J).AND.XN.GE.XLIMS1(1,J).AND.
     .        YN.LE.YLIMS2(1,J).AND.YN.GE.YLIMS1(1,J).AND.
     .        ZN.LE.ZLIMS2(1,J).AND.ZN.GE.ZLIMS1(1,J)
        ELSE
          XMS1=XN*PS13(1,J)+YN*PS13(2,J)+ZN*PS13(3,J)+P1A(J)
          XLS1=XN*PS23(1,J)+YN*PS23(2,J)+ZN*PS23(3,J)+P2A(J)
          LGJ=XMS1.GE.0..AND.XLS1.GE.0..AND.XMS1+XLS1.LE.1.
          IF (RLB(J).GE.4.AND..NOT.LGJ) THEN
            XMS2=XN*PS24(1,J)+YN*PS24(2,J)+ZN*PS24(3,J)+P1B(J)
            XLS2=XN*PS34(1,J)+YN*PS34(2,J)+ZN*PS34(3,J)+P2B(J)
            LGJ=XMS2.GE.0..AND.XLS2.GE.0..AND.XMS2+XLS2.LE.1.
            IF (RLB(J).GE.5.AND..NOT.LGJ) THEN
              XMS3=XN*PS35(1,J)+YN*PS35(2,J)+ZN*PS35(3,J)+P1C(J)
              XLS3=XN*PS45(1,J)+YN*PS45(2,J)+ZN*PS45(3,J)+P2C(J)
              LGJ=XMS3.GE.0..AND.XLS3.GE.0..AND.XMS3+XLS3.LE.1.
            ENDIF
          ENDIF
        ENDIF
        IF (RLBNOT(J)) LGJ=.NOT.LGJ
        IF (NLPRCS(J).AND.LGJ) LCNDEXP=.TRUE.
        IF (LTSTCXP) GOTO 100
        TMX=TMA
        IF (LGJ) GOTO 70
        GOTO 100
C
C   A1*TMX+A2=0
C
40      TMA=-A2/A1
        IF (NLTRC) WRITE (iunout,*) 
     .    'TIMEA AT 40, J,TMA,A1,A2 ',J,TMA,A1,A2
!PB        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=-A2
        IF (RLB(J)) 26,16,36
C
C   A2*TMX+A3=0
C
50      IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
        IF (NLTRC) WRITE (iunout,*) 
     .    'TIMEA AT 50, J,TMA,A2,A3 ',J,TMA,A2,A3
!PB        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2
        IF (RLB(J)) 26,16,36
C
C  A1LM(J).NE.0
C
60      CONTINUE
        A2=A1LM(J)*VX+A2LM(J)*VY+A3LM(J)*VZ
        A3=A0LM(J)+A1LM(J)*X+A2LM(J)*Y+A3LM(J)*Z
        IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
!PB        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
        XN=ALM(J)+BLM(J)*YN+CLM(J)*ZN
        TUP=TMA
        ICOUNT=2
        IF (RLB(J)) 28,18,38
C
C  A2LM(J).NE.0
C
63      CONTINUE
        A2=A1LM(J)*VX+A2LM(J)*VY+A3LM(J)*VZ
        A3=A0LM(J)+A1LM(J)*X+A2LM(J)*Y+A3LM(J)*Z
        IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
!PB        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2
        XN=X+TMA*VX
        ZN=Z+TMA*VZ
        YN=ALM(J)+BLM(J)*XN+CLM(J)*ZN
        TUP=TMA
        ICOUNT=2
        IF (RLB(J)) 28,18,38
C
C  A3LM(J).NE.0
C
66      CONTINUE
        A2=A1LM(J)*VX+A2LM(J)*VY+A3LM(J)*VZ
        A3=A0LM(J)+A1LM(J)*X+A2LM(J)*Y+A3LM(J)*Z
        IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
!PB        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2
        XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=ALM(J)+BLM(J)*XN+CLM(J)*YN
        TUP=TMA
        ICOUNT=2
        IF (RLB(J)) 28,18,38
C
C
C   DATA RETURNED TO CALLING PROGRAM
C   TENTATIVELY FOR SURFACE NO. J
C
70      TL=TMX+TADD
        TMIN=TMX
        XR=XN
        YR=YN
        ZR=ZN
        NNR=NN
        VXJ=VX
        VZJ=VZ
        MASURF=J
        SG=SIGN(1._DP,WR)
        IF (NLTRC) THEN
          WRITE (iunout,*) 'TIMEA, TL,XR,YR,ZR,NNR,MASURF,SG '
          WRITE (iunout,*)         TL,XR,YR,ZR,NNR,MASURF,SG
        ENDIF
C
C  LOOP OVER SURFACE-INDEX FINISHED
C
100   CONTINUE
C
C **********************************************************************
C
C  IF NO INTERSECTION FOUND, RETURN
      IF (MASURF.EQ.0) RETURN
C  INTERSECTION AT SURFACE NO. MASURF
C  IF NOT TRANSPARENT, RETURN
      IF (ILIIN(MASURF).GT.0) RETURN
C  IF TRANSPARENT BUT WRONG SIDE, RETURN
      IF (ILSIDE(MASURF)*SG.LT.0) RETURN
C
      IF (NLTRC) WRITE (iunout,*) 'NOT RETURNED FROM TIMEA, OTHER LOOP '
C  ILIIN=0, CONTINUE WITH ANOTHER LOOP IN SUBR. TIMEA
      IF (ILIIN(MASURF).EQ.0) THEN
        X=XR
        XS=X
        Y=YR
        YS=Y
        Z=ZR
        ZS=Z
C
        TADD=TL
        NLLLI=MASURF
        NN=NNR
        NNTCLS=NN
        VX=VXJ
        VXS=VXJ
        VZ=VZJ
        VZS=VZJ
        GOTO 1000
      ELSEIF (ILIIN(MASURF).LT.0) THEN
C  TRANSPARENT, BUT SWITCH AND/OR SURFACE TALLIES
        RETURN
      ENDIF
C
      END
C ===== SOURCE: timep.f
C
C
      SUBROUTINE TIMEP (ZRAD)
C
C   CALCULATE TIME SEGMENTS IN Y- OR POLOIDAL MESH
C
C   INPUT:
C
C       X00,Y00,Z00: STARTING POINT FOR THIS TRACK
C       ZRAD   =
C       NRCELL = RADIAL CELL NUMBER, IN WHICH THIS TRACK OF LENGTH ZRAD
C                                    IS PERFORMED
C              = 0, IF TRACK OUTSIDE RADIAL MESH
C                   IF MRSURF .NE. 0 THEN
C                   REENTRY FOUND AT RADIAL SURFACE MRSURF.
C                   OTHERWISE: REENTRY AT NON DEFAULT
C                   POLOIDAL SURFACES IS SEARCHED FOR IN THIS CALL
C       NTCELL =
C       ITCELL =
C       NPCELL = POLOIDAL CELL INDEX OF POINT X00,Y00,Z00
C                AT WHICH THIS TRACK OF LENGTH ZRAD STARTS
C  WARNING: IF NLSRFY, NPCELL MAY BE WRONG
C       IPCELL =
C
C       NCOUT,BLPD,NCOUNT
C
C   OUTPUT:
C
C       X00,Y00,Z00: END POINT FOR THIS TRACK
C       NCOU   = TOTAL NUMBER OF CELLS, WHICH THE CURRENT TRACK OF
C                LENGTH ZRAD HAS CROSSED
C       CLPD(I)= LENGTH OF THE PART NO. I (I=1,NCOU)
C       JUPC(I)= CELL NUMBER IN P-GRID FOR TRACK I
C       LUPC(I)= SURFACE NUMBER IN P-GRID AN THE END OF TRACK I
C       MUPC(I)= ORIENTATION OF TRACK I
C       NPCELL = LAST POLOIDAL CELL INDEX OF X00,Y00,Z00 POINT,
C                ON WHICH THIS TRACK OF LENGTH ZRAD ENDS
C       IPCELL = NUMBER OF FINAL POLOIDAL CELL NPCELL
C                IF TRACK ORIGINATED INSIDE STANDARD MESH
C              = NUMBER OF POLOIDAL CELL, AT WHICH REENTRY WAS FOUND
C                ON RADIAL OR TOROIDAL SURFACE
C                IF TRACK ORIGINATED OUTSIDE  STANDARD MESH
C       IRCELL = NUMBER OF RADIAL CELL, AT WHICH REENTRY AT POLOIDAL
C                SURFACE MPSURF WAS FOUND. OTHERWISE: UNMODIFIED.
C
C              = 0  IF TRACK COMPLETELY OUTSIDE STANDARD MESH
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE COMPRT
      USE COMSOU
      USE CLGIN

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: ZRAD
      REAL(DP) :: HELP, GS, GC, F, XNEN, DUM, SUM, BB, DY, X0TEST,
     .          Y0TEST, V1, V2, T1, ZRADS, ZRD, WIN1, X000, Z0T,
     .          TWIN1, XT, Y000, TSS, X0T, Y0T
      INTEGER :: ILLZ, IANP, I, IENP, JHELP, ISW, JJC, J1, NRCLLP,
     .           LHELP, J2, INCY, JN, ICOU, NYSAVE, NHELP,
     .           LEARCA, I1, MPTEST, IR, IRSAVE, NCOUPE, ISTS, IADD,
     .           ICOUT, NCPEN, IPOLGS, NCPAN, J, LEARC2, NJC, LEARC1,
     .           IST, ITEST, JSH, NN ,IN
c slmod begin - gfortran
      INTEGER, ALLOCATABLE :: NCOUNS(:)
      LOGICAL, ALLOCATABLE :: LCUTY(:),LCUTX(:)
      LOGICAL :: lnincz
c
c      INTEGER :: NCOUNS(N2ND+N3RD)
c      LOGICAL :: LCUTY(N2NDPLG), LCUTX(N1ST), lnincz
c slmod end
      SAVE
c slmod begin - gfortran
      IF (.NOT.ALLOCATED(NCOUNS)) ALLOCATE(NCOUNS(N2ND+N3RD))
      IF (.NOT.ALLOCATED(LCUTY)) ALLOCATE(LCUTY(N2NDPLG))
      IF (.NOT.ALLOCATED(LCUTX)) ALLOCATE(LCUTX(N1ST))
c slmod end
C     
      IF (NLTRC) THEN
        CALL LEER(1)
        IF (NRCELL.GT.0) THEN
          WRITE (iunout,*) 'TIMEP FROM INSIDE, NPANU ', NPANU
          WRITE (iunout,*) 'ZRAD,NRCELL,NPCELL '
          WRITE (iunout,*)  ZRAD,NRCELL,NPCELL
        ELSE
          WRITE (iunout,*) 'TIMEP FROM OUTSIDE, NPANU ', NPANU
          WRITE (iunout,*) 'MRSURF,MTSURF,ZRAD '
          WRITE (iunout,*)  MRSURF,MTSURF,ZRAD
        ENDIF
      ENDIF
C
      ICOUT=1
      IADD=0
      ZRADS=ZRAD
      IPOLGS=IPOLGN
      ncpan=0
      ncpen=0
C
      IF (NLTOR) THEN
C       NCOUT=NCOUT
        ZRAD=BLPD(1)
        NTCELL=NCOUNT(1)
        IF (NCOUT.GT.1) IPOLGN=0
        IF (NLTRC) WRITE (iunout,*) 
     .    'WG. NLTOR: ZRAD,NTCELL ',ZRAD,NTCELL
      ELSE
        NCOUT=1
C       ZRAD=ZRAD
C       NTCELL=1
      ENDIF
C
10000 CONTINUE
C
CDR NOV.99: ADDED BECAUSE OF FOLION OPTION, WITH DEFAULT B-FIELD (IN Z DIRECTION
C
      IF (ABS(VELZ).EQ.1.D0) THEN
        NCOUP=1
        ALPD(NCOUP)=ZRAD
        JUPC(NCOUP)=NPCELL
C       X00=X00+ZRAD*VELX
C       Y00=Y00+ZRAD*VELY
        IPCELL=NPCELL
        MPSURF=0
        GOTO 5000
      ENDIF
C
      IF (LEVGEO.EQ.3) THEN
C
        IF (NRCELL.GT.0) GOTO 10
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT POLOIDAL SURFACES
C
        NCOUP=0
        ZRD=ZRAD
        DO 3 ISTS=1,NSTSI
          MPTEST=INUMP(ISTS,2)
          IF (MPTEST.NE.0) THEN
C  TEST POLOIDAL SURFACE NO. MPTEST FOR REENTRY
            DO 2 IR=1,NR1STM
              I1=IR+1
              V1=(YPOL(IR,MPTEST)-Y00)*VELX-(XPOL(IR,MPTEST)-X00)*VELY
              V2=(YPOL(I1,MPTEST)-Y00)*VELX-(XPOL(I1,MPTEST)-X00)*VELY
              LCUTX(IR)=V1*V2.LE.0.
2           CONTINUE
            DO 4 IR=1,NR1STM
              IF (LCUTX(IR)) THEN
                T1=((XPOL(IR,MPTEST)-X00)*VVTY(IR,MPTEST)-
     .              (YPOL(IR,MPTEST)-Y00)*VVTX(IR,MPTEST))
     .             /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)
                IF (NLTRC) WRITE (iunout,*) 'IR,MPTEST,T1 ',IR,MPTEST,T1
                IF (T1.LT.0..OR.T1.GE.ZRD) GOTO 4
                IF (NLTRC) WRITE (iunout,*) 'VALID INTERSECTION AT T1= '
     .                                      ,T1
                NCOUP=1
                LUPC(NCOUP)=MPTEST
                IRSAVE=IR
                MUPC(NCOUP)=SIGN(1._DP,VELX*PPLNX(IR,MPTEST)+
     .                                 VELY*PPLNY(IR,MPTEST))
                JUPC(NCOUP)=1
                ZRD=T1
              ENDIF
4           CONTINUE
          ENDIF
3       CONTINUE
        IF (NCOUP.GT.0.AND.LUPC(MAX(1,NCOUP)).NE.MPSURF) THEN
C  REENTRY FOUND, REDUCE ZRAD TO T1
C  NCOUP=1 AT THIS POINT
          NCOUPE=1
          MPSURF=LUPC(NCOUPE)
          IPOLGN=LUPC(NCOUPE)
          NINCY=MUPC(NCOUPE)
          IRCELL=IRSAVE
          ZRAD=ZRD
          ISRFCL=0
          ALPD(NCOUPE)=ZRAD
          MRSURF=0
          MTSURF=0
          MASURF=0
          NINCX=0
          NINCZ=0
          IF (NLTRC) THEN
            WRITE (iunout,*) 'REENTRY FOUND, MPSURF,ZRAD = ',MPSURF,ZRAD
            WRITE (iunout,*) 'IRCELL ',IRCELL
          ENDIF
          GOTO 31
        ELSE
C  NO REENTRY FOUND
          NCOUP=1
          JUPC(1)=1
          ALPD(1)=ZRAD
          IPCELL=IPOLGN
          MPSURF=0
          IF (NLTRC) THEN
            WRITE (iunout,*) 'NO REENTRY INTO POLOIDAL GRID FOUND '
          ENDIF
          X00=X00+ZRAD*VELX
          Y00=Y00+ZRAD*VELY
          GOTO 5000
        ENDIF
C
10      CONTINUE
C
C  PARTICLE INSIDE STANDARD MESH, RADIAL CELL NO. NRCELL
C
        IF (NCOUP.EQ.0) THEN
          WRITE (iunout,*) 'ERROR IN TIMEP: NCOUP=0'
          RETURN
        ENDIF
        IF (NLTRC) THEN
          WRITE (iunout,*) ' TIMEP IN NEIGHBOR PART '
          WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
          WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
          WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
          WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
        ENDIF
        NLSRFY=.FALSE.
        IF ((icout == 1) .and.
     .    (SQRT((X0-X00)**2+(Y0-Y00)**2).GT.EPS10)) THEN
          DO J=1,NCOUP
            ALPD(J)=ALPD(J)-ZT
          ENDDO
        ENDIF
C
C  ACCOUNT FOR OTHER SURFACES INSIDE POLOIDAL MESH
C
        lnincz=.false.
        if (ityp==3) lnincz=(zrad < alpd(1))
        IF (ALPD(NCOUP).GT.ZRAD) THEN
          DO J=1,NCOUP
            IF (ALPD(J).GT.ZRAD) THEN
              ncpan=j+1
              ncpen=ncoup+1
              do jsh=ncoup,j+1,-1
                alpd(jsh+1)=alpd(jsh)-zrad
                jupc(jsh+1)=jupc(jsh)
                lupc(jsh+1)=lupc(jsh)
                mupc(jsh+1)=mupc(jsh)
              end do
              alpd(ncpan)=alpd(j)-zrad
              jupc(ncpan)=jupc(j)
              lupc(ncpan)=lupc(j)
              mupc(ncpan)=mupc(j)
              ALPD(J)=ZRAD
              IPOLGN=JUPC(J)
              NCOUP=J
              LUPC(NCOUP)=0
              MUPC(NCOUP)=0
              GOTO 4711
            ENDIF
          ENDDO
4711      CONTINUE
        ENDIF
C
C   ADJUST ALPD AND ACCOUNT FOR "NONDEFAULT" POLOIDAL SURFACES
C
        TSS=0.
        IST=0
        DO J=1,NCOUP
          ALPD(J)=ALPD(J)-TSS
          TSS=TSS+ALPD(J)
C
          ITEST=INMP2I(NRCELL,LUPC(J),0)
          IN=ITEST+NLIM
          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
C
C  TRACK ENDS ON ONE OF THE NON DEFAULT POLOIDAL SURFACES
C
            IF (NLTRC) THEN
              WRITE (iunout,*) ' TRACK TERMINATED'
              WRITE (iunout,*) ' ITEST,ILIIN ',ITEST,ILIIN(IN)
            ENDIF
            NCOUPE=J
            MPSURF=LUPC(NCOUPE)
            IPOLGN=LUPC(NCOUPE)
            NINCY=MUPC(NCOUPE)
            ZRAD=TSS
            ISRFCL=0
            MASURF=0
            nincx=0
            nincz=0
            mrsurf=0
            mtsurf=0
            NCOUP=NCOUPE
            GOTO 311

          ELSEIF (ityp==3) THEN
C
C  STOP TRACK ANYHOW
C
            IF (NLTRC) THEN
              WRITE (iunout,*) ' TRACK TERMINATED'
              WRITE (iunout,*) ' ITYP,ILIIN ',ITYP,ILIIN(IN)
            ENDIF
            NCOUPE=J
            MPSURF=LUPC(NCOUPE)
            IF (LUPC(NCOUPE) /= 0) IPOLGN=LUPC(NCOUPE)
            NINCY=MUPC(NCOUPE)
            ZRAD=TSS
            if (1.-zrad/(zrads+eps60) > eps10) then
              ISRFCL=0
              MASURF=0
            end if
!pb
            if (lnincz) then         ! toroidal surface
              nincx=0
              nincy=0
              mrsurf=0
              mpsurf=0
            elseif ((ityp==3).and.(nincx.ne.0)) then   ! radial surface
              nincy=0
              nincz=0
              mpsurf=0
              mtsurf=0
            else                     ! poloidal surface
              nincx=0
              nincz=0
              mrsurf=0
              mtsurf=0
            end if
!pb
            NCOUP=NCOUPE
            GOTO 311
          ENDIF
C
        ENDDO
C
C  LAST CELL, TRACK DOES NOT END ON A POLOIDAL SURFACE
C
C  INDEX IPOLGN OF LAST CELL NOT KNOWN ?
        IF (IPOLGN.EQ.0) THEN
          X0T=X00+VELX*ZRAD
          Y0T=Y00+VELY*ZRAD
          NN=LEARC1(X0T,Y0T,Z0T,IPOLGN,
     .              NRCELL,NRCELL,.FALSE.,.FALSE.,
     .              NPANU,'TIMEP       ')
        ENDIF
C
        MPSURF=0
        NINCY=0
C
311     CONTINUE
        IF (NLTRC) THEN
          WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
          WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
          WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
          WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
        ENDIF
C
        X00=X00+ZRAD*VELX
        Y00=Y00+ZRAD*VELY
        NPCELL=JUPC(NCOUP)
        IPCELL=NPCELL
        GOTO 5000
C
C
30      CONTINUE
C
C  LAST CELL, TRACK DOES NOT END ON A POLOIDAL SURFACE
C
C  INDEX OF LAST CELL NOT KNOWN (E.G. DUE TO ADD. SURFACE) ?
CDR  ERROR: FALLS NLSRFY, GGFLS NPCELL FALSCH
        IF (IPOLGN.EQ.0) THEN
CDR       IF (NCOUP.EQ.0) THEN
CDR         IPOLGN=NPCELL
CDR       ELSE
            X0T=X00+VELX*ZRAD
            Y0T=Y00+VELY*ZRAD
            IPOLGN=LEARC2(X0T,Y0T,NRCELL,NPANU,'TIMEP       ')
CDR       ENDIF
        ENDIF
C
        NCOUPE=NCOUP+1
        JUPC(NCOUPE)=IPOLGN
        LUPC(NCOUPE)=0
        MUPC(NCOUPE)=0
        ALPD(NCOUPE)=ZRAD-TSS
        MPSURF=0
        NINCY=0
C
31      CONTINUE
        NCOUP=NCOUPE
        IF (NLTRC) THEN
          WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
          WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
          WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
          WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
        ENDIF
C
        X00=X00+ZRAD*VELX
        Y00=Y00+ZRAD*VELY
        NPCELL=JUPC(NCOUP)
        IPCELL=NPCELL
        GOTO 5000
C
      ELSEIF (LEVGEO.EQ.2) THEN
        IF (NLCRC) THEN
C
C
C  DISTANCE TO NEXT X- OR RADIAL SURFACE KNOWN?
C
C       IF (TS.GE.1.D30) GOTO 1000
C
C  YES! ZRAD IS THE DISTANCE TRAVELLED IN X- OR RADIAL CELL NO. NRCELL
C
          NCOUP=1
          IF (NRCELL.LT.1.OR.NRCELL.GT.NR1STM) THEN
            X00=X00+VELX*ZRAD
            Y00=Y00+VELY*ZRAD
            WIN1=MOD(ATAN2(Y00,X00)+PI2A-PSURF(1),PI2A)+PSURF(1)
            NPCELL=WIN1/YDF*DBLE(NP2NDM)+1.
            GOTO 5000
          ENDIF
C
C  THE OLD POLOIDAL CELL INDEX IS: NPCELL
C  FIND THE NEW CELL INDEX : NJC
C
          X000=X00+VELX*ZRAD
          Y000=Y00+VELY*ZRAD
          WIN1=MOD(ATAN2(Y000,X000)+PI2A-PSURF(1),PI2A)+PSURF(1)
          NJC=WIN1/YDF*DBLE(NP2NDM)+1.
C
          TWIN1=0.
          NCOUP=0
          IF (NJC.EQ.NPCELL) GOTO 150
C
C   FIND ORIENTATION IN THETA-GRID
C
          XT=-Y00*VELX+X00*VELY
          NINCY=1
          IF (XT.LT.0.) NINCY=-1
C
C   CONTRIBUTION TO EACH THETA-CELL
C   NPCELL : STARTINDEX
C   JJC    : SURFACEINDEX
C   J1     : CELLINDEX
C   NJC    : ENDINDEX
          J1=NPCELL
100       JJC=J1
          IF (NINCY.EQ.1) JJC=JJC+1
C   TIMESTEP FROM X00,Y00 TO THETA-SURFACE, THETA=WIN
          GS=SINPH(JJC)
          GC=COSPH(JJC)
          XNEN=VELX*GS-VELY*GC
          F=(Y00*GC-X00*GS)/(XNEN+EPS60)
          NCOUP=NCOUP+1
          JUPC(NCOUP)=J1
          LUPC(NCOUP)=JJC
          ALPD(NCOUP)=F-TWIN1
          TWIN1=F
C
          J1=J1+NINCY
          IF (J1.EQ.0) J1=NP2NDM
          IF (J1.EQ.NP2ND) J1=1
!pb          IF (J1.NE.NJC) GOTO 100
          IF ((ityp.ne.3).and.(J1.NE.NJC)) GOTO 100
C
C   LAST THETA-CELL
C
150       NCOUP=NCOUP+1
          JUPC(NCOUP)=NJC
          ALPD(NCOUP)=ZRAD-TWIN1
          X00=X000
          Y00=Y000
          NPCELL=NJC
          MPSURF=0
          NINCY=0
          GOTO 5000
C
        ELSE
C
C  NEW PART: NOT NLCRC
C  PARTICLE INSIDE STANDARD MESH
C
          NCOUP=0
          NRCLLP=NRCELL+1
C
C   SEARCH FOR ALL POSSIBLE INTERSECTIONS WITHIN THE RADIAL CELL NRCELL
C
          DO 111 I=1,NP2ND
111         LCUTY(I)=.FALSE.
C
          DO 112 J=1,NP2ND
            V1=(YPOL(NRCELL,J)-Y00)*VELX-(XPOL(NRCELL,J)-X00)*VELY
            V2=(YPOL(NRCLLP,J)-Y00)*VELX-(XPOL(NRCLLP,J)-X00)*VELY
            LCUTY(J)=V1*V2.LE.0.
            IF (NLTRC) THEN
              IF (LCUTY(J)) WRITE (iunout,*) 'LCUTY=TRUE FOR ',J
            ENDIF
112       CONTINUE
          IF (NLSRFY) THEN
            LCUTY(MPSURF)=.FALSE.
            NLSRFY=.FALSE.
C  PSURF(1)=PSURF(NP2ND)
            IF (MPSURF.EQ.1) LCUTY(NP2ND)=.FALSE.
            IF (MPSURF.EQ.NP2ND) LCUTY(1)=.FALSE.
          ENDIF
C
          IANP=ILLZ(NP2ND,LCUTY,1)+1
          IENP=NP2ND-ILLZ(NP2ND,LCUTY,-1)
C
C   COMPUTE THE FLIGHT TIMES TO THE INTERSECTION POINTS
C
        DO 114 I=IANP,IENP
          IF (LCUTY(I)) THEN
            T1=((XPOL(NRCELL,I)-X00)*VVTY(NRCELL,I)-
     .          (YPOL(NRCELL,I)-Y00)*VVTX(NRCELL,I))
     .         /(VELX*VVTY(NRCELL,I)-VELY*VVTX(NRCELL,I)+EPS60)
            IF (NLTRC) WRITE (iunout,*) 'I,T1 ',I,T1
            IF (T1.LT.0..OR.T1.GE.ZRAD) GOTO 114
            IF (NLTRC) WRITE (iunout,*) 'VALID INTERSECTION AT T1= ',T1
            NCOUP=NCOUP+1
            LUPC(NCOUP)=I
            MUPC(NCOUP)=SIGN(1._DP,VELX*PPLNX(NRCELL,I)+
     .                             VELY*PPLNY(NRCELL,I))
            IF (MUPC(NCOUP).EQ.-1) THEN
              JUPC(NCOUP)=I
              ALPD(NCOUP)=T1
              IF (I.EQ.NP2ND) NCOUP=NCOUP-1
            ELSEIF (MUPC(NCOUP).EQ.1) THEN
              JUPC(NCOUP)=I-1
              ALPD(NCOUP)=T1
              IF (I.EQ.1) NCOUP=NCOUP-1
            ENDIF
          ENDIF
114     CONTINUE
C
C   REARRANGE THE FLIGHT TIMES IN ASCENDING ORDER
C
115     ISW=0
        DO 120 J=1,NCOUP-1
          IF (ALPD(J).GT.ALPD(J+1)) THEN
            ISW=ISW+1
            HELP=ALPD(J)
            ALPD(J)=ALPD(J+1)
            ALPD(J+1)=HELP
            JHELP=JUPC(J)
            JUPC(J)=JUPC(J+1)
            JUPC(J+1)=JHELP
            LHELP=LUPC(J)
            LUPC(J)=LUPC(J+1)
            LUPC(J+1)=LHELP
            NHELP=MUPC(J)
            MUPC(J)=MUPC(J+1)
            MUPC(J+1)=NHELP
          ENDIF
120     CONTINUE
        IF (ISW.GT.0.AND.NCOUP.GT.2) GOTO 115
        IF (NLTRC.AND.NCOUP.GT.0) THEN
          WRITE (iunout,*) ' NACH SORTIEREN '
          WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
          WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
          WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
          WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
        ENDIF
C
        DO 125 J=1,NCOUP-1
          IF (ABS(ALPD(J+1)-ALPD(J)).LE.EPS30) THEN
            IF (JUPC(J).LE.0.OR.JUPC(J).GE.NP2ND) THEN
              IF (NLTRC) THEN
                WRITE (iunout,*) ' VERTAUSCHE ALPD(',J,') UND (',J+1,')'
                WRITE (iunout,*) ' ALPD = ',ALPD(J),ALPD(J+1)
                WRITE (iunout,*) ' JUPC = ',JUPC(J),JUPC(J+1)
              ENDIF
              HELP=ALPD(J)
              ALPD(J)=ALPD(J+1)
              ALPD(J+1)=HELP
              JHELP=JUPC(J)
              JUPC(J)=JUPC(J+1)
              JUPC(J+1)=JHELP
              LHELP=LUPC(J)
              LUPC(J)=LUPC(J+1)
              LUPC(J+1)=LHELP
              NHELP=MUPC(J)
              MUPC(J)=MUPC(J+1)
              MUPC(J+1)=NHELP
            ENDIF
          ENDIF
125     CONTINUE
C
C   ADJUST ALPD AND ACCOUNT FOR "NONDEFAULT" POLOIDAL SURFACES
C
        TSS=0.
        IST=0
        DO 130 J=1,NCOUP
          ALPD(J)=ALPD(J)-TSS
          TSS=TSS+ALPD(J)
C
          ITEST=INMP2I(NRCELL,LUPC(J),0)
          IN=ITEST+NLIM
!pb          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
          IF ((ityp==3).or.(ITEST.NE.0.AND.ILIIN(IN).NE.0)) THEN
C
C  TRACK ENDS ON ONE OF THE NON DEFAULT POLOIDAL SURFACES
C
            IF (NLTRC) THEN
              WRITE (iunout,*) ' TRACK TERMINATED'
              WRITE (iunout,*) ' ITEST,ILIIN ',ITEST,ILIIN(IN)
            ENDIF
            NCOUPE=J
            MPSURF=LUPC(NCOUPE)
            IPOLGN=LUPC(NCOUPE)
            NINCY=MUPC(NCOUPE)
            ZRAD=TSS
            ISRFCL=0
            NINCX=0
            NINCZ=0
            MRSURF=0
            MTSURF=0
            MASURF=0
            GOTO 131
          ENDIF
C
130     CONTINUE
C
C  LAST CELL, TRACK DOES NOT END ON A POLOIDAL SURFACE
C
C  INDEX OF LAST CELL NOT KNOWN (E.G. DUE TO ADD. SURFACE) ?
CDR  ERROR: FALLS NLSRFY, GGFLS NPCELL FALSCH
        IF (IPOLGN.EQ.0) THEN
CDR       IF (NCOUP.EQ.0) THEN
CDR         IPOLGN=NPCELL
CDR       ELSE
            X0T=X00+VELX*ZRAD
            Y0T=Y00+VELY*ZRAD
            IPOLGN=LEARC2(X0T,Y0T,NRCELL,NPANU,'TIMEP       ')
CDR       ENDIF
        ENDIF
C
        NCOUPE=NCOUP+1
        JUPC(NCOUPE)=IPOLGN
        LUPC(NCOUPE)=0
        MUPC(NCOUPE)=0
        ALPD(NCOUPE)=ZRAD-TSS
        MPSURF=0
        NINCY=0
C
131       CONTINUE
          NCOUP=NCOUPE
          IF (NLTRC) THEN
            WRITE (iunout,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
            WRITE (iunout,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
            WRITE (iunout,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
            WRITE (iunout,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
          ENDIF
C
          X00=X00+ZRAD*VELX
          Y00=Y00+ZRAD*VELY
          NPCELL=JUPC(NCOUP)
          IPCELL=NPCELL
          GOTO 5000
C
        ENDIF
C
      ELSEIF (LEVGEO.EQ.1) THEN
C
C  IDENTICAL, UP TO NAMES, TO TOROIDAL GRID PART, NLTRZ BLOCK
C
        IF (NRCELL.GT.0) GOTO 2900
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT POLOIDAL SURFACES
C
        NCOUP=0
        ZRD=ZRAD
        BB=VELY+EPS60
        NYSAVE=1
        IF (VELY.LT.0.D0) NYSAVE=-1
        DO 2903 ISTS=1,NSTSI
          MPTEST=INUMP(ISTS,2)
          IF (MPTEST.NE.0) THEN
C  TEST POLOIDAL SURFACE NO. MPTEST FOR REENTRY
C  TIME FROM Y00 TO PSURF
            DY=PSURF(MPTEST)-Y00
            F=DY/BB
            IF (NLTRC) WRITE (iunout,*) 'MPTEST,F,DY ',MPTEST,F,DY
            IF (F.LE.ZRD.AND.F.GT.0.D0) THEN
              X0TEST=X00+VELX*F
              IF (X0TEST.GE.RSURF(1).AND.X0TEST.LE.RSURF(NR1ST)) THEN
                IRSAVE=LEARCA(X0TEST,RSURF,1,NR1ST,1,'TIMEP 1    ')
                NCOUP=1
                JUPC(NCOUP)=1
                LUPC(NCOUP)=MPTEST
                MUPC(NCOUP)=NYSAVE
                ZRD=F
              ENDIF
            ENDIF
          ENDIF
2903    CONTINUE
        IF (NCOUP.GT.0.AND.LUPC(MAX(1,NCOUP)).NE.MPSURF) THEN
C  REENTRY FOUND, REDUCE ZRAD TO F
C  NCOUP=1 AT THIS POINT
          NCOUP=1
          MPSURF=LUPC(NCOUP)
          NINCY=MUPC(NCOUP)
          IRCELL=IRSAVE
          ZRAD=ZRD
          ISRFCL=0
          ALPD(NCOUP)=ZRAD
          MRSURF=0
          MTSURF=0
          MASURF=0
          NINCX=0
          NINCZ=0
          IF (NLTRC) THEN
            WRITE (iunout,*) 'REENTRY FOUND, MPSURF,ZRAD = ',MPSURF,ZRAD
            WRITE (iunout,*) 'IRCELL ',IRCELL
          ENDIF
          Y00=Y00+VELY*ZRD
          NPCELL=JUPC(NCOUP)
          IPCELL=NPCELL
          GOTO 5000
        ELSE
C  NO REENTRY FOUND
          NCOUP=1
          JUPC(1)=1
          ALPD(1)=ZRAD
          IF (MRSURF.GT.0) THEN
C  CHECK VALID RANGE ON MRSURF
            Y0TEST=Y00+VELY*ZRAD
            IF (NLTRC) WRITE (iunout,*) 'CHECK VALID RANGE: Y0TEST ',
     .        Y0TEST
            IF (Y0TEST.GE.PSURF(1).AND.Y0TEST.LE.PSURF(NP2ND)) THEN
              IPCELL=LEARCA(Y0TEST,PSURF,1,NP2ND,1,'TIMEP 2     ')
            ELSE
              MRSURF=0
              MTSURF=0
              NINCX=0
              NINCZ=0
            ENDIF
          ENDIF
          MPSURF=0
          NINCY=0
          IF (NLTRC) THEN
            WRITE (iunout,*) 'NO REENTRY FOUND '
          ENDIF
          Y00=Y00+ZRD*VELY
          GOTO 5000
        ENDIF
C
C  PARTICLE IN STANDARD MESH, RADIAL CELL NRCELL
C
2900    CONTINUE
        Y000=Y00+VELY*ZRAD
        IF (NLTRC) WRITE (iunout,*) 'Y000,Y00,ZRAD ',Y000,Y00,ZRAD
C
        DUM=0.
        NCOUP=1
C
C  J1: CELL INDEX
C  J2: SURFACE INDEX
        J1=NPCELL
        IF (VELY.LT.0.) THEN
          INCY=0
          NINCY=-1
          IF (NLSRFY) J1=MPSURF-1
        ELSE
          INCY=1
          NINCY=1
          IF (NLSRFY) J1=MPSURF
        ENDIF
        J2=J1+INCY
C
        NLSRFY=.FALSE.
C
3000    CONTINUE
        IF (J2.LE.0.OR.J2.GT.NP2ND) THEN
          WRITE (iunout,*) 'ERROR IN TIMEP ',J2,J1,VELY
          CALL EXIT_OWN(1)
        ENDIF
C  TIME FROM Y00 TO PSURF
        IF (MPSURF.EQ.J2) THEN
          J1=J1+NINCY
          J2=J1+INCY
          GOTO 3000
        ENDIF
        DY=(PSURF(J2)-Y00)
!pb reduce zrad if the trajectory comes too close to the next standard surface
!pb without hitting it
        IF (ABS(ABS(DY/ZRAD)-1._DP) <= EPS10) THEN
          ZRAD=(1._DP-EPS6)*ZRAD
          ZRADS=ZRAD
        END IF
        BB=VELY+EPS60
        F=DY/BB
        IF (NLTRC) WRITE (iunout,*) 'J2,F,DY ',J2,F,DY
        IF (F.LE.ZRAD) THEN
          JUPC(NCOUP)=J1
          LUPC(NCOUP)=J2
          ALPD(NCOUP)=F-DUM
          DUM=F
C  STOP HISTORY AT NON DEFAULT STANDARD SURFACE J2
          ITEST=INMP2I(0,J2,0)
          IN=ITEST+NLIM
!pb          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
          IF ((ityp==3).or.(ITEST.NE.0.AND.ILIIN(IN).NE.0)) THEN
            ZRAD=F
            ISRFCL=0
            NPCELL=J1
            MPSURF=J2
            MRSURF=0
            NINCX=0
            NINCZ=0
            MTSURF=0
            MASURF=0
            IPCELL=NPCELL
            Y00=PSURF(J2)
            GOTO 5000
          ENDIF
C  NEXT CELL
          J1=J1+NINCY
          J2=J1+INCY
          NCOUP=NCOUP+1
          GOTO 3000
        ENDIF
C
C  LAST CELL
C
3100    CONTINUE
        IF (NLTRC) WRITE (iunout,*) 'LAST CELL ',ZRAD,DUM
        IF (MPSURF.EQ.0) THEN
          NPCELL=J1
        ELSE
          Y0T=Y00+VELY*ZRAD
          NPCELL=LEARCA(Y0T,PSURF,1,NP2ND,1,'TIMEP 3     ')
        ENDIF
        MPSURF=0
        JUPC(NCOUP)=J1
        ALPD(NCOUP)=ZRAD-DUM
        IPCELL=NPCELL
        Y00=Y000
        GOTO 5000
C
      ENDIF
C
5000  CONTINUE
      DO 5100 J=1,NCOUP
        CLPD(IADD+J)=ALPD(J)
        NUPC(IADD+J)=(JUPC(J)-1)+(NTCELL-1)*NP2T3
        NCOUNP(IADD+J)=JUPC(J)
        NCOUNS(IADD+J)=NTCELL
        IF (ALPD(J).LE.0..OR.JUPC(J).LE.0.OR.JUPC(J).GE.NP2ND) THEN
          WRITE (iunout,*) 'ERROR DETECTED IN TIMEP '
          WRITE (iunout,*) 'J,IADD+J,ALPD,JUPC ',
     .                      J,IADD+J,ALPD(J),JUPC(J)
        ENDIF
        IF (NLTRC) WRITE (iunout,*) 'TIMEP ',
     .      J+IADD,CLPD(J+IADD),NUPC(J+IADD),NCOUNP(J+IADD)
5100  CONTINUE
C
      IF (ICOUT.LT.NCOUT.AND.MPSURF.EQ.0) THEN
        ICOUT=ICOUT+1
        IPOLGN=0
        IF (ICOUT.EQ.NCOUT) IPOLGN=IPOLGS
        ZRAD=BLPD(ICOUT)
        if (ncpan.gt.0) then
          jn=0
          do j=ncpan,ncpen
            jn=jn+1
            alpd(jn)=alpd(j)
            jupc(jn)=jupc(j)
            lupc(jn)=lupc(j)
            mupc(jn)=mupc(j)
          end do
        endif
        NTCELL=NCOUNT(ICOUT)
        IADD=IADD+NCOUP
        ncoup=jn
        IF (NLTRC) WRITE (iunout,*) 
     .    'NEXT TOR. CELL: ZRAD,NTCELL,IADD ',ZRAD,NTCELL,IADD
        GOTO 10000
      ENDIF
C
      NCOU=IADD+NCOUP
C
      SUM=0.
      DO 5110 ICOU=1,NCOU
        SUM=SUM+CLPD(ICOU)
        NCOUNT(ICOU)=NCOUNS(ICOU)
C       WRITE (iunout,*) 'ICOU,NCOUNT,CLPD ',
C    .                    ICOU,NCOUNT(ICOU),CLPD(ICOU)
5110  CONTINUE
      IF (MPSURF.EQ.0.AND.ABS(SUM-ZRADS).GT.EPS10) THEN
        WRITE (iunout,*) 'ERROR IN TIMEP: NPANU,SUM,ZRADS ',
     .                               NPANU,SUM,ZRADS
        RETURN  ! to avoid job crash in long runs
!        CALL EXIT_OWN(1)
      ENDIF
C
      ZRAD=SUM
      RETURN
C
991   CONTINUE
      WRITE (iunout,*) 
     .  'REENTRANCE FROM VACUUM REGION AND LEVGEO=1 IN TIMEP'
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: timer.f
C
C  FULL EIRENE GEOMETRY BLOCK  (GEO3D)
C
C
      SUBROUTINE TIMER (PT)
C
C  THIS SUBROUTINE CALCULATES INTERSECTION TIMES IN THE STANDARD
C  MESH "RSURF" (X- OR RADIAL DIRECTION)
C
C  INPUT:
C       NRCELL = CELL NUMBER FOR WHICH NEXT INTERSECTION
C                TIME IS TO BE CALCULATED (I.E. NOT NECESSARLY THE
C                CELL WHICH CONTAINS THE STARTING POINT X0,Y0,Z0)
C                NRCELL=0 IF PARTICLE OUTSIDE STANDARD MESH
C       NJUMP = 0 MEANS: THIS IS THE FIRST CALL OF TIMER FOR THIS TRACK
C                        IN THIS CASE , NLSRFX MUST BE KNOWN
C           NLSRFX = .TRUE. :PARTICLE ON A SURFACE, IN THIS CASE
C                            THE SUBROUTINE NEEDS 'MRSURF'
C                            MRSURF = NUMBER OF THIS SURFACE
C
C           NLSRFX = .FALSE.:PARTICLE NOT ON A SURFACE
C
C           X0,Y0,Z0 = STARTING POINT OF THIS TRACK
C           VELX,VELY,VELZ = VELOCITY OF PARTICLE
C       NJUMP = 1  X0,Y0,Z0,VELX,VELY,VELZ ARE THE SAME
C                  AS IN THE PREVIOUS CALL
C       NJUMP = 2  ONLY VELX,VELY,VELZ ARE THE SAME
C                  AS IN THE PREVIOUS CALL, I.E. PARTICLE HAS
C                  BEEN MOVED BUT VELOCITY HAS NOT BEEN CHANGED
C                  (TO BE WRITTEN)
C  IF (LEVGEO=1 OR LEVGEO=2):
C       TIMINT NE 0. MEANS: INTERSECTION TIME IS KNOWN FROM AN EARLIER
C                CALL AND ITS VALUE IS TIMINT.
C                OTHERWISE (TIMINT=0) IT HAS TO BE CALCULATED IN THIS
C                CALL
C  IF (LEVGEO=3):
C       TIMINT NE 0. MEANS: INTERSECTION TIME IS KNOWN FROM AN EARLIER
C                CALL AND IS TO BE FOUND IN THE ARRAYS TIMPOL,IIMPOL.
C                OTHERWISE (TIMINT=0) IT HAS TO BE CALCULATED IN THIS
C                CALL
C  IF (LEVGEO=4):
C       TIMINT NE 0. MEANS: NOTHING
C  IF (LEVGEO=5):
C       TIMINT NE 0. MEANS: NOTHING
C  IF (LEVGEO=6):
C       TIMINT NE 0. MEANS: NOTHING
C
C  OUTPUT :
C       NJUMP = 1
C       MRSURF = INDEX OF NEXT SURFACE ALONG TRACK
C              = 0 IF NO NEXT SURFACE IS FOUND
C       PT = TIME TO REACH THIS SURFACE
C          = 1.D30 IF NO NEXT SURFACE IS FOUND
C       TIMINT(MRSURF) NE 0 INDICATES FURTHER INTERSECTION TIMES FOUND
C                      IN THIS CALL, WHICH MAY BE USED IN A LATER CALL
C       NINCX = INDEX FOR DIRECTION IN GRID "RSURF": +1 OR -1
C             = 0 IF NO NEXT SURFACE FOUND
C       NRCELL:  NUMBER OF FINAL RADIAL CELL (IE. NOT MODIFIED)
C       IRCELL: NOT NEEDED ON RADIAL SURFACE. SET IRCELL=NRCELL
C
C  ADDITIONALLY IF NLPLG:
C  INPUT  :
C       IPOLG  = INDEX ON POLYGON OF INITIAL POINT X=(X0,Y0,Z0)
C  OUTPUT :
C       IPOLGN = INDEX ON POLYGON OF THE POINT X+PT*VEL
C         ( IF NO VALID INTERSECTION FOUND IN THIS CALL, THEN
C           IPOLGN=IPOLGO IS RETURNED, THE INDEX OF THE LAST
C           VALID POINT OF INTERSECTION FOUND IN EARLIER CALLS
C           (WHICH MAY BE THE INPUT VALUE IPOLG ITSELF) )
C
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CTSURF
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTETRA
      USE COMPRT
      USE CLGIN
      USE CTRIG

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: PT
      REAL(DP) :: A(3,3), B(3), AB(3,3)
      REAL(DP) :: SIG1, SIG2, DET, XMU, XETA, SARRUS, AX, AY, V, TIMTET,
     .          ZZ, RICHTY, HELP, TM, RICHTX, YY, XX, S3, TIMT, TIM,
     .          S2, A1, A2, A3, B1, B2, B3, S1, V2, V3, ZR, ZEP1, 
     .          ZSQRT, PS, VELYQ, XX0, VVELX, VELXQ, ZT1, ZT2, YVY,
     .          ZC1, DXA, TST, XA, ZB2, ZAB, ZAB2, ZB, Z0TEST, 
     .          Y0TEST, ZA, T, PT1, PT2, PT3, PT4, V1, ESURF, 
     .          XTEST, YTEST, X0SURF, DSRF, Y0Q, T1, T2, T3, T4, PNORMI
      INTEGER :: IRICH(2,4), ITSIDE(3,4)
      INTEGER :: IZELLO, NTIMT, IPOLGOO, IOB, I1, I2, ISW, KAN, KEN,
     .           ILLZ, IHELP, IZELL, NTMS, MXSF, NTMZ, NRMSRF, ICOS,
     .           IERR, IRS, NEWCEL, ITET, IT, IL, IS, NRI, MS, IR,
     .           ICALL, ITFRST, ISTS, MMSURF, ICOUP, J, K, I, JPOL,
     .           MPOL, IPOLGO, LEARC2, IP, ICELLR, MSAVE, ITRI, 
     .           ISTS_CELL
c slmod begin - gfortran
      INTEGER, ALLOCATABLE :: ITRINO(:), ISIDNO(:)
      INTEGER :: NSTS_CELL
      LOGICAL, ALLOCATABLE :: LCUT(:)
c
c      INTEGER, ALLOCATABLE, SAVE :: ITRINO(:), ISIDNO(:)
c      INTEGER, SAVE :: NSTS_CELL
c      LOGICAL :: LCUT(N2NDPLG)
c slmod end
      LOGICAL :: LNGB1, LNGB2, LNGB3, LNGB4,
     .           LCT1, LCT2, LCT3, LCT4, BITGET

      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/
      DATA IRICH / 1, -3,
     .             4,  1,
     .             5,  2,
     .             6,  3 /
      DATA ICALL /0/
      DATA ITFRST /0/
      SAVE
c slmod begin - gfortran
      IF (.NOT.ALLOCATED(LCUT)) ALLOCATE(LCUT(N2NDPLG))
c slmod end
C     
      IF (NLTRC) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'TIMER, INIT.: NRCELL,NLSRFX,MRSURF,TT,TL'
        WRITE (iunout,*)                NRCELL,NLSRFX,MRSURF,TT,TL
        WRITE (iunout,*) '              NPCELL,NLSRFY,MPSURF'
        WRITE (iunout,*)                NPCELL,NLSRFY,MPSURF
      ENDIF
C
      IRCELL=NRCELL
      PT=1.D30
      IF ((LEVGEO <= 4) .AND. (ABS(VELZ).EQ.1.D0)) THEN
        IPOLGN=IPOLG
        NCOUP=1
        ALPD(NCOUP)=PT
        LUPC(NCOUP)=0
        MUPC(NCOUP)=0
        JUPC(NCOUP)=IPOLGN
        RETURN
      ENDIF
C
C-------------------------------------------------------------------
      IF (LEVGEO.GT.1) GOTO 100
C
C****SLAB-MODEL IN X DIRECTION
C
      IF (ABS(VELY).EQ.1.D0) THEN
        IPOLGN=IPOLG
        RETURN
      ENDIF
C
      IF (NJUMP.EQ.0) THEN
        XA=X0
        IF (NLSRFX) XA=RSURF(MRSURF)
        NINCX=1
        IF (VELX.LT.0.) NINCX=-1
C
C   RESET BRANCH MARKER
        NJUMP=1
        NLSRFX=.FALSE.
      ENDIF
C
      MRSURF=0
      IF (NRCELL.GT.0) THEN
        IR=NRCELL
        IF (NINCX.EQ.1) IR=IR+1
        DXA=RSURF(IR)-XA
        TST=DXA/(VELX+EPS60)
        IF (NLTOR.AND.NLTRZ) THEN
          Z0TEST=Z0+TST*VELZ
          IF (ZSURF(1).GT.Z0TEST.OR.ZSURF(NT3RD).LT.Z0TEST) GOTO 5
        ENDIF
        IF (NLPOL) THEN
          Y0TEST=Y0+TST*VELY
          IF (PSURF(1).GT.Y0TEST.OR.PSURF(NP2ND).LT.Y0TEST) GOTO 5
        ENDIF
        PT=TST
        MRSURF=IR
5       CONTINUE
      ELSE
C  TRY TO FIND REENTRY SURFACE. CHECK ONLY NON DEFAULT RADIAL SURFACES
        DO 10 IR=1,NR1ST
          IF (INMP1I(IR,0,0).NE.0) THEN
            DXA=RSURF(IR)-XA
            TST=DXA/(VELX+EPS60)
            IF (TST.LE.0..OR.TST.GT.PT) GOTO 10
            IF (NLTOR.AND.NLTRZ) THEN
              Z0TEST=Z0+TST*VELZ
              IF (ZSURF(1).GT.Z0TEST.OR.ZSURF(NT3RD).LT.Z0TEST) GOTO 10
            ENDIF
            IF (NLPOL) THEN
              Y0TEST=Y0+TST*VELY
              IF (PSURF(1).GT.Y0TEST.OR.PSURF(NP2ND).LT.Y0TEST) GOTO 10
            ENDIF
            PT=TST
            MRSURF=IR
          ENDIF
10      CONTINUE
        IF (MRSURF.EQ.0) NINCX=0
      ENDIF
      IF (NLTRC) THEN
        WRITE (iunout,*) 'TIMER, OUT: PT,MRSURF,NINCX'
        WRITE (iunout,*)              PT,MRSURF,NINCX
      ENDIF
      RETURN
C
C
100   CONTINUE
C
C---------------------------------------------------------------------
      IF (LEVGEO.GT.2) GOTO 6000
C
      IF (NLELL) GOTO 1000
C
C****CIRCULAR MESH   (CONCENTRIC)
C
      IF (NLTRC) THEN
        WRITE (iunout,*) 'NJUMP,NINCX,NRCELL '
        WRITE (iunout,*)  NJUMP,NINCX,NRCELL
      ENDIF
      IF (NJUMP.EQ.0) THEN
        ZA=VELX*VELX+VELY*VELY
        ZB=X0*VELX+Y0*VELY
        ZB2=ZB*ZB
        ZAB=-ZB/(ZA+EPS60)
        ZAB2=ZA/(ZB2+EPS60)
C
C  TEST FOR DIRECTION
        NINCX=1
        IF (ZB.LT.0) NINCX=-NINCX
C
        IF (NLSRFX) THEN
          ZC1=RQ(MRSURF)
C  IN THIS CASE: DON'T TRUST THE VALUE OF NRCELL. RECOMPUTE FROM MRSURF
          IF (NINCX.EQ.1) NRCELL=MRSURF
          IF (NINCX.EQ.-1) NRCELL=MRSURF-1
        ELSE
          ZC1=X0*X0+Y0*Y0
        ENDIF
      ENDIF
C
      IF (NRCELL.GT.0) THEN
C  PARTICLE INSIDE STANDARD MESH
C  FIND NEXT SURFACE MRSURF
        NRI=0
        MRSURF=NRCELL
        IF (NINCX.EQ.1) MRSURF=MRSURF+1
        GOTO 204
      ENDIF
C
C  PARTICLE OUTSIDE STANDARD MESH
C  FIND NEXT SURFACE MRSURF
      NRI=1
      PS=PT
      MS=0
      IS=0
200   DO 205 IR=NRI,NR1ST
        IF (INMP1I(IR,0,0).NE.0) THEN
          MRSURF=IR
          GOTO 204
        ENDIF
205   CONTINUE
      PT=PS
      MRSURF=MS
      RETURN
C
C  CHECK SURFACE: MRSURF
C
204   IF (MRSURF.EQ.1) GO TO 201
      IF (TIMINT(MRSURF).GT.0.0) GO TO 303
203   ZR=RQ(MRSURF)
C  CHECK FOR ROOT
      ZEP1=ZAB2*(ZC1-ZR)
      IF (ZEP1.LT.1.0) GO TO 202
C  NO ROOT - PATH IN OTHER DIRECTION. THIS MUST BE THE
C  OUTWARD DIRECTION NOW, DUE TO CONVEXITY OF THE MESH
      IF (NRI.GT.0) GOTO 310
201   CONTINUE
      NINCX=1
      MRSURF=MRSURF+1
      IF (NJUMP.EQ.0) GO TO 203
      GO TO 303
202   CONTINUE
C
C  RESET BRANCH MARKER
      NJUMP=1
C
C  INTERSECTION TIMES
      ZSQRT=1.0+SQRT(1.0-ZEP1)
      ZT1=ZAB*ZEP1/ZSQRT
      ZT2=ZAB*ZSQRT
      NLSRFX=.FALSE.
C
CL              3.         RETURN RESULT
C
300   CONTINUE
C
C  CHECK SIGN OF RESULT
      IF(ZEP1.GT.0.0) GO TO 301
C  ONE ROOT NEGATIVE OR ZERO
      PT=MAX(ZT1,ZT2)
      GOTO 310
C
301   CONTINUE
      PT=ZT1
      TIMINT(MRSURF)=ZT2
      NIMINT = NIMINT+1
      IIMINT(NIMINT) = MRSURF
      GOTO 310
C
C  ROOT ALREADY KNOWN
303   PT=TIMINT(MRSURF)
      GOTO 310
C
310   CONTINUE
      IF (NRI.EQ.0) THEN
        RETURN
      ELSE
        IF (PT.GT.0.AND.PT.LT.PS) THEN
          IS=MRSURF-MAX(0,NINCX)
          MS=MRSURF
          PS=PT
        ENDIF
        NRI=IR+1
        GOTO 200
      ENDIF
C
1000  CONTINUE
C
C****ELLIPTICAL MESH  (NOT NECESSARILY CONCENTRIC OR CONFOCAL)
C
      IF (NRCELL.LE.0) THEN
        WRITE (iunout,*) 'NRCELL.LE.0 IN TIMER: NOT READY '
        CALL EXIT_OWN(1)
      ENDIF

C  PARTICLE OUTSIDE GRID OPTION NOT READY, SEE NLCRC
      NRI=0

C  COEFFICIENTS OF QUADRATIC
      IF (NJUMP.EQ.0) THEN
        YVY=Y0*VELY
        VELXQ=VELX*VELX
        VELYQ=VELY*VELY
        XX0=X0
        VVELX=VELX
        IF (NLSRFX) THEN
          X0SURF=X0-EP1(MRSURF)
          DSRF=ELLQ(MRSURF)
          ZR=RQ(MRSURF)
          Y0Q=(ZR-X0SURF*X0SURF)*DSRF
          MSAVE=MRSURF
C  IN THIS CASE: DON'T TRUST THE VALUE OF NRCELL. RECOMPUTE FROM MRSURF
          ZB=X0SURF*VELX+Y0*VELY/ELLQ(MRSURF)
C  TEST FOR DIRECTION. CASE: NLSRFX, MRSURF KNOWN
          NINCX=1
          IF (ZB.LT.0)     NINCX =-NINCX
          IF (NINCX.EQ.1)  NRCELL=MRSURF
          IF (NINCX.EQ.-1) NRCELL=MRSURF-1
C  NEXT SURFACE
          MRSURF=NRCELL
          IF (NINCX.EQ.1) MRSURF=MRSURF+1
          GOTO 1100
        ELSE
          MSAVE=0
          Y0Q=Y0*Y0
C   TEST FOR DIRECTION. CASE:.NOT.NLSRFX, NRCELL KNOWN
          NINCX=-1
          MRSURF=NRCELL
          ZB=(XX0-EP1(MRSURF))*VVELX+YVY/ELLQ(MRSURF)
          IF(ZB.LE.0.) GOTO 1200
          NINCX=1
          MRSURF=MRSURF+1
          GOTO 1100
        ENDIF
      ENDIF
C
2000  CONTINUE
      MRSURF=NRCELL
      IF (NINCX.EQ.1) MRSURF=MRSURF+1
1100  ZB=(XX0-EP1(MRSURF))*VVELX+YVY/ELLQ(MRSURF)
C
1200  IF(MRSURF.EQ.1) GO TO 2010
      IF(TIMINT(MRSURF).GT.0.0) GO TO 3030
      ESURF=XX0-EP1(MRSURF)
      DSRF=ELLQ(MRSURF)
      ZA=VELXQ+VELYQ/DSRF
      ZB2=ZB*ZB
      ZAB=-ZB/(ZA+EPS60)
      ZAB2=ZA/(ZB2+EPS60)
      ZC1=ESURF*ESURF+Y0Q/DSRF
      ZR=RQ(MRSURF)
      ZEP1=ZAB2*(ZC1-ZR)
C  CHECK FOR ROOT
      IF(ZEP1.LT.1.0) GO TO 2020
C  NO ROOT - PATH IN OTHER DIRECTION
C  MUST BE OUTWARD, BECAUSE OF CONVEXITY OF MESH
2010  CONTINUE
      NINCX=1
      MRSURF=MRSURF+1
      IF (NJUMP.EQ.0.) GO TO 1100
      GO TO 3030
2020  CONTINUE
C
C  RESET BRANCH MARKER
      NJUMP=1
C
C  INTERSECTION TIMES
      ZSQRT=1.0+SQRT(1.0-ZEP1)
      IF (MRSURF.EQ.MSAVE) THEN
        ZT1=0.
        ZEP1=0.
      ELSE
        ZT1=ZAB*ZEP1/ZSQRT
      ENDIF
      ZT2=ZAB*ZSQRT
      NLSRFX=.FALSE.
C
C---------------------------------------------------------------------
CL              3.         RETURN RESULT
C
3000  CONTINUE
C
C  CHECK SIGN OF RESULT
      IF(ZEP1.GT.0.0) GO TO 3010
C  ONE ROOT NEGATIVE OR ZERO
      PT=MAX(ZT1,ZT2)
      GOTO 3100
C
3010  CONTINUE
C  BOTH ROOTS POSITIVE - RETURN ONE AND SAVE OTHER
      PT=ZT1
      TIMINT(MRSURF)=ZT2
      NIMINT = NIMINT+1
      IIMINT(NIMINT) = MRSURF
      GOTO 3100
C
C  ROOT ALREADY KNOWN
3030  PT=TIMINT(MRSURF)
      GOTO 3100

C
3100  CONTINUE
      XTEST=X0+PT*VELX
      YTEST=Y0+PT*VELY
      IF (NLPOL) IPOLGN=LEARC2(XTEST,YTEST,NRCELL,NPANU,'TIMER')
      IF (NLTRC) 
     .  WRITE (iunout,*) 'IPOLGN,NRCELL FROM TIMER ',IPOLGN,NRCELL
      IF (NRI.EQ.0) THEN
        RETURN
      ELSE
        WRITE (iunout,*) ' OPTION NRI > 0 NOT READY in TIMER '
        CALL EXIT_OWN(1)
CPB        IF (PT.GT.0.AND.PT.LT.PS) THEN
CPB          IS=MRSURF-MAX(0,NINCX)
CPB          MS=MRSURF
CPB          PS=PT
CPB        ENDIF
CPB        NRI=IR+1
CPB        GOTO 2000
      ENDIF
C
C     POLYGONS
C
6000  CONTINUE
C
      IF (LEVGEO.GT.3) GOTO 8000
C
C  POLYGON MESH
C
      IF (NRCELL.NE.0) THEN
        IF (NLTRC) WRITE (iunout,*) ' TIMER: IN NEIGHBOR PART'
        LNGB1=.TRUE.
        LNGB2=.TRUE.
        LNGB3=.TRUE.
        LNGB4=.TRUE.
C  SET INDICES OF STARTING CELL
        IR=NRCELL
        IP=NPCELL
        IF (NJUMP.EQ.1) THEN
          IR=ICELLR
          IP=IPOLGN
          IF (NINCX.EQ.1) LNGB1=.FALSE.
          IF (NINCX.EQ.-1) LNGB3=.FALSE.
        ENDIF
        IF (NLSRFX) THEN
          NLSRFX=.FALSE.
          IF (NRCELL.EQ.MRSURF) THEN
            LNGB1=.FALSE.
          ELSE
            LNGB3=.FALSE.
          ENDIF
        ENDIF
        IF (NLSRFY) THEN
          IF (IPOLG.EQ.MPSURF) THEN
            IF (NPCELL.EQ.MPSURF) THEN
              LNGB4=.FALSE.
              IP=NGHPLS(2,IR,MPSURF)
            ELSE
              LNGB2=.FALSE.
              IP=NGHPLS(4,IR,MPSURF)
            ENDIF
          ELSE
            LNGB2=.FALSE.
          ENDIF
        ENDIF
C
        NCOUP=0
C  CALCULATE INTERSECTIONS OF FLIGHT WITH CELL BOUNDARIES
6001    CONTINUE
        IF (NLTRC) WRITE (iunout,*) ' IR,IP,LNGB1,LNGB2,LNGB3,LNGB4',
     .                           IR,IP,LNGB1,LNGB2,LNGB3,LNGB4
        T1=-1.D30
        T2=-1.D30
        T3=-1.D30
        T4=-1.D30
        PT1=-1.D30
        PT2=-1.D30
        PT3=-1.D30
        PT4=-1.D30
        IF (LNGB1)
     .  T1=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .     (VELX*VPLY(IR,IP)-VELY*VPLX(IR,IP)+EPS60)
        IF (LNGB2)
     .  T2=((XPOL(IR,IP+1)-X0)*VELY-(YPOL(IR,IP+1)-Y0)*VELX)/
     .      (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
        IF (LNGB3)
     .  T3=((XPOL(IR+1,IP)-X0)*VELY-(YPOL(IR+1,IP)-Y0)*VELX)/
     .     (VELX*VPLY(IR+1,IP)-VELY*VPLX(IR+1,IP)+EPS60)
        IF (LNGB4)
     .  T4=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .      (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)
        IF (NLTRC) WRITE (iunout,*) ' T1,T2,T3,T4 ',T1,T2,T3,T4
        LCT1 = T1.GE.0.D0 .AND. T1.LE.1.D0
        LCT2 = T2.GE.0.D0 .AND. T2.LE.1.D0
        LCT3 = T3.GE.0.D0 .AND. T3.LE.1.D0
        LCT4 = T4.GE.0.D0 .AND. T4.LE.1.D0
C  CALCULATE TIME OF FLIGHT FROM STARTING POINT X=(X0,Y0,Z0)
C  TO THE BOUNDARY OF THE ACTUELL CELL
        IF (ABS(VELX).GT.ABS(VELY)) THEN
          IF (LCT1) PT1=(XPOL(IR,IP)-X0+VPLX(IR,IP)*T1)/VELX
          IF (LCT2) PT2=(XPOL(IR,IP+1)-X0+VVTX(IR,IP+1)*T2)/VELX
          IF (LCT3) PT3=(XPOL(IR+1,IP)-X0+VPLX(IR+1,IP)*T3)/VELX
          IF (LCT4) PT4=(XPOL(IR,IP)-X0+VVTX(IR,IP)*T4)/VELX
        ELSE
          IF (LCT1) PT1=(YPOL(IR,IP)-Y0+VPLY(IR,IP)*T1)/VELY
          IF (LCT2) PT2=(YPOL(IR,IP+1)-Y0+VVTY(IR,IP+1)*T2)/VELY
          IF (LCT3) PT3=(YPOL(IR+1,IP)-Y0+VPLY(IR+1,IP)*T3)/VELY
          IF (LCT4) PT4=(YPOL(IR,IP)-Y0+VVTY(IR,IP)*T4)/VELY
        ENDIF
        LCT1 = LCT1 .AND. PT1.GE.0.D0
        LCT2 = LCT2 .AND. PT2.GE.0.D0
        LCT3 = LCT3 .AND. PT3.GE.0.D0
        LCT4 = LCT4 .AND. PT4.GE.0.D0
        IF (NLTRC) WRITE (iunout,*) ' LCT1,LCT2,LCT3,LCT4 ',
     .                           LCT1,LCT2,LCT3,LCT4
        IF (NLTRC) WRITE (iunout,*) ' PT1,PT2,PT3,PT4 ',PT1,PT2,PT3,PT4
        T=MAX(PT1,PT2,PT3,PT4)
        IF (NLTRC) WRITE (iunout,*) ' T = ',T
C  IF INTERSECTION WITH POLOIDAL BOUNDARY CONTINUE WITH NEIGHBORING CELL
        IF (LCT2.OR.LCT4) THEN
          LNGB1=.TRUE.
          LNGB2=.TRUE.
          LNGB3=.TRUE.
          LNGB4=.TRUE.
          NCOUP=NCOUP+1
          ALPD(NCOUP)=T
          IF (LCT2) THEN
            LUPC(NCOUP)=IP+1
            MUPC(NCOUP)=1
            JUPC(NCOUP)=IP
            IP=NGHPOL(2,IR,IP)
            LNGB4=.FALSE.
          ENDIF
          IF (LCT4) THEN
            LUPC(NCOUP)=IP
            MUPC(NCOUP)=-1
            JUPC(NCOUP)=IP
            IP=NGHPOL(4,IR,IP)
            LNGB2=.FALSE.
          ENDIF
          ISTS=INMP2I(IR,LUPC(NCOUP),0)
!pb          IF ((.not.NLPOL.or.ISTS.eq.0).and.ip.ne.0) goto 6001
          IF (ityp.ne.3.and.(.not.NLPOL.or.ISTS.eq.0).and.ip.ne.0)
     .       goto 6001
C  NO NEIGHBORING CELL: PARTICLE HAS HIT A POLOIDAL BOUNDARY OF THE MESH
          MRSURF=0
          PT=1.D30
          NINCX=0
        ELSEIF (LCT1.OR.LCT3) THEN
C  INTERSECTION WITH RADIAL CELL BOUNDARY FOUND
          PT=T
          IPOLGN=IP
          LNGB1=.TRUE.
          LNGB2=.TRUE.
          LNGB3=.TRUE.
          LNGB4=.TRUE.
          NCOUP=NCOUP+1
          ALPD(NCOUP)=T
          LUPC(NCOUP)=0
          MUPC(NCOUP)=0
          JUPC(NCOUP)=IP
          IF (LCT1) THEN
            MRSURF=IR
            NINCX=-1
            ICELLR=IR-1
          ELSE
            MRSURF=IR+1
            ICELLR=IR+1
            NINCX=1
          ENDIF
        ELSE
C  NO INTERSECTION FOUND
          IF (NLSRFY.AND.NJUMP.EQ.0) THEN
C  PLAY SAVE: TRY ONCE AGAIN, IF PARTICLE ON POL. SURFACE
            WRITE (iunout,*) 
     .        ' NO INTERSECTION IN TIMER. TRY ONCE AGAIN '
            MMSURF=MSURF
            IF (MSURF.GT.NLIM) MMSURF=-MSURF+NLIM
            WRITE (iunout,*) 'NPANU, MSURF ',NPANU,MMSURF
            IF (.NOT.LNGB4) THEN
              IP=NGHPLS(4,IR,MPSURF)
              LNGB4=.TRUE.
              LNGB2=.FALSE.
            ELSEIF (.NOT.LNGB2) THEN
              IP=NGHPLS(2,IR,MPSURF)
              LNGB2=.TRUE.
              LNGB4=.FALSE.
            ENDIF
            LNGB1=.TRUE.
            LNGB3=.TRUE.
            NJUMP=1
            GOTO 6001
          ENDIF
          WRITE (iunout,*) ' ERROR: NO INTERSECTION FOUND IN TIMER '
          WRITE (iunout,*) ' NPANU: ',NPANU
          MRSURF=0
          PT=1.D30
          NINCX=0
          ICELLR=0
          IPOLGN=IP
        ENDIF
        NJUMP=1
        IF (NLTRC) THEN
          WRITE (iunout,*) ' PT,MRSURF,NINCX,IPOLGN ',
     .                  PT,MRSURF,NINCX,IPOLGN
          WRITE (iunout,*) 'NCOUP ',NCOUP
          DO ICOUP=1,NCOUP
            WRITE (iunout,*) 'ICOUP,ALPD(ICOUP) ',ICOUP,ALPD(ICOUP)
          ENDDO
        ENDIF
        RETURN
      ENDIF
C
C  PARTICLE OUTSIDE STANDARD MESH, IN ADDITIONAL CELL NACELL
C  NRCELL=0
C
      IF (NLSRFX) THEN
        NLSRFX=.FALSE.
        JPOL=IPOLG
        MPOL=MRSURF
      ELSE
        JPOL=0
        MPOL=0
      ENDIF
C
      IF (NJUMP.EQ.0) IPOLGO=IPOLG
C
      PT=1.D30
      MRSURF=0
      IPOLGN=IPOLGO
C
      DO 6100 I=1,NR1ST
        ISTS=INMP1I(I,0,0)
C  SURFACE I IS NOT A NON DEFAULT RADIAL STANDARD SURFACE
        IF (ISTS.EQ.0) GOTO 6100
        IF (TIMINT(I).EQ.0) THEN
          TIMINT(I)=1
          NIMINT = NIMINT+1
          IIMINT(NIMINT) = I
          NTIM(I)=0
C
C   SEARCH FOR ALL POSSIBLE INTERSECTIONS WITHIN THE CELL
C
          DO 6011 J=1,NRPLG
6011        LCUT(J)=.FALSE.
C
          DO 6012 J=1,NPPLG
            DO 6012 K=NPOINT(1,J),NPOINT(2,J)-1
              V1=(YPOL(I,K+1)-Y0)*VELX-(XPOL(I,K+1)-X0)*VELY
              V2=(YPOL(I,K)-Y0)*VELX-(XPOL(I,K)-X0)*VELY
              LCUT(K)=V1*V2.LE.0.
6012      CONTINUE
C
          IF (I.EQ.MPOL) LCUT(JPOL)=.FALSE.
          KAN=ILLZ(NRPLG,LCUT,1)+1
          KEN=NRPLG-ILLZ(NRPLG,LCUT,-1)
C
C   COMPUTE THE FLIGHT TIMES TO THE INTERSECTION POINTS
C
          DO 6013 K=KAN,KEN
            IF (LCUT(K)) THEN
              T1=((XPOL(I,K)-X0)*VPLY(I,K)-(YPOL(I,K)-Y0)*VPLX(I,K))
     .           /(VELX*VPLY(I,K)-VELY*VPLX(I,K)+EPS60)
              IF (T1.GT.0.) THEN
                NTIM(I)=NTIM(I)+1
                TIMPOL(I,NTIM(I))=T1
                IIMPOL(I,NTIM(I))=K
              ENDIF
            ENDIF
6013      CONTINUE
C
6015      ISW=0
          DO 6020 J=1,NTIM(I)-1
            IF (TIMPOL(I,J).LT.TIMPOL(I,J+1)) THEN
              ISW=ISW+1
              HELP=TIMPOL(I,J)
              TIMPOL(I,J)=TIMPOL(I,J+1)
              TIMPOL(I,J+1)=HELP
              IHELP=IIMPOL(I,J)
              IIMPOL(I,J)=IIMPOL(I,J+1)
              IIMPOL(I,J+1)=IHELP
            ENDIF
6020      CONTINUE
          IF (ISW.GT.0.AND.NTIM(I).GT.2) GOTO 6015
          IF (NLTRC) THEN
            WRITE (iunout,*) ' SURFACE NO. = ',I
            WRITE (iunout,*) ' TIMPOL ',(TIMPOL(I,J),J=1,NTIM(I))
            WRITE (iunout,*) ' IIMPOL ',(IIMPOL(I,J),J=1,NTIM(I))
          ENDIF
C
        ENDIF
C
C  FIND PT, MRSURF, IPOLGN
        IF (NTIM(I).GT.0) THEN
          IF (TIMPOL(I,NTIM(I)).LT.PT) THEN
            MRSURF=I
            PT=TIMPOL(I,NTIM(I))
            IPOLGN=IIMPOL(I,NTIM(I))
          ENDIF
        ENDIF
6100  CONTINUE
C
C  FIND NINCX
      NINCX=0
      IF (MRSURF.NE.0) THEN
        NTIM(MRSURF)=NTIM(MRSURF)-1
        NINCX=SIGN(1._DP,VELX*PLNX(MRSURF,IPOLGN)+
     .        VELY*PLNY(MRSURF,IPOLGN))
      ENDIF
C
      NJUMP=1
      IPOLGO=IPOLGN
C
      IF (NLTRC) WRITE (iunout,*) 'PT,MRSURF,IPOLGN,NINCX ',
     .                        PT,MRSURF,IPOLGN,NINCX
      RETURN
C
C
8000  CONTINUE
C
      IF (LEVGEO.GT.4) GOTO 10000
C
C   FINITE ELEMENT DISCRETISATION
C
      IF (NLTRC) THEN
        WRITE (iunout,*) ' TIMER,NJUMP,NRCELL ',NJUMP,NRCELL
        WRITE (iunout,*) '       NLSRFX,IPOLG ',NLSRFX,IPOLG
      ENDIF
C
      IF (NJUMP.EQ.0) THEN
8001    CONTINUE
        IZELL = NRCELL
        IPOLGO=0
        IF (NLSRFX) THEN
          IPOLGO=IPOLG
        ENDIF
C     ELSEIF (NJUMP.EQ.1) THEN
      ENDIF

      IF (NRCELL > 0) THEN
        IF (NLTRC) THEN
          WRITE (iunout,*) ' IZELL,IPOLGO ',IZELL,IPOLGO
          WRITE (iunout,*) ' X0,Y0 ',X0,Y0
          WRITE (iunout,*) ' VELX,VELY ',VELX,VELY
          CALL LEER(1)
        ENDIF
        XX = X0
        YY = Y0
        TM=0.
        IF (IZELL.EQ.0) GOTO 9999
C
8020    CONTINUE
        DO 8010,J=1,3
          IF (J.NE.IPOLGO) THEN
            RICHTX = VTRIX(J,IZELL)
            RICHTY = VTRIY(J,IZELL)
            AX = XTRIAN(NECKE(J,IZELL))-XX
            AY = YTRIAN(NECKE(J,IZELL))-YY
            V  = (AX*VELY-AY*VELX)/(RICHTY*VELX-RICHTX*VELY+EPS60)
            IF (V.GE.0..AND.V.LE.1.) THEN
              IF (ABS(VELX).GT.ABS(VELY)) THEN
                T=(AX+V*RICHTX)/VELX
              ELSE
                T=(AY+V*RICHTY)/VELY
              ENDIF
              IF (T .GT. 0.) THEN
                XX = XX + T * VELX
                YY = YY + T * VELY
                TM = TM + T
                TIMINT(IZELL) = TM
                NIMINT = NIMINT+1
                IIMINT(NIMINT) = IZELL
C  NUMMER DER GESCHNITTENEN SEITE DES ALTEN DREIECKS
                NTIM(IZELL) = J
                IF (NLTRC) THEN
                  WRITE (iunout,*) 'IZELL,XX,YY,J ',IZELL,XX,YY,J
                ENDIF
C  ZELLENNUMMER DES NEUEN DREIECKS
                IZELLO=IZELL
                IZELL = NCHBAR(J,IZELLO)
                IF (IZELL .EQ. 0) GOTO 8050
C  SEITENNUMMER DES NEUEN DREIECKS
                IPOLGO = NSEITE(J,IZELLO)
                IF (NLTRC) THEN
                  WRITE(iunout,*) 'GEHE IN ZELLE ',IZELL,' TM= ',TM
                  WRITE(iunout,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
                ENDIF
                GOTO 8050
              ENDIF
            ENDIF
          ENDIF
8010    CONTINUE
C
        IF (NLTRC) 
     .    WRITE (iunout,*) ' NO INTERSECTION FOUND IN TRIANGLE '
        PT=1.D30
        NINCX=0
C  DURCHSUCHE NACHBARDREIECK
        IF (NLSRFX) THEN
          IZELL=NCHBAR(IPOLG,NRCELL)
c slmod begin 
          IF (IZELL.EQ.0) THEN
            WRITE(6,*) 'ERROR: IZELL=0, PARTICLE ABANDONED'
            WRITE(6,*) NPANU,NRCELL
            WRITE(0,*) 'ERROR: IZELL=0, PARTICLE ABANDONED'
            WRITE(0,*) NPANU,NRCELL
            RETURN
          ENDIF
c slmod end
          IPOLGO=NSEITE(IPOLG,NRCELL)
          NRCELL=IZELL
          NLSRFX=.FALSE.
          GOTO 8020
        ENDIF
C  KEIN SCHNITTPUNKT GEFUNDEN UND NLSRFX=.FALSE.
        RETURN
C
C
C   THIS IS DONE FOR NJUMP=0 AND NJUMP=1
C
8050  CONTINUE
      NLSRFX=.FALSE.
      NJUMP=1
      MRSURF=NRCELL
      PT=TIMINT(MRSURF)
      TIMINT(MRSURF)=0.
      NIMINT = NIMINT+1
      IIMINT(NIMINT) = MRSURF
      IPOLGN=NTIM(NRCELL)
      ISTS=ABS(INMTI(IPOLGN,NRCELL))
      IF (ISTS.EQ.0) THEN
        NINCX=NCHBAR(IPOLGN,NRCELL)-NRCELL
      ELSEIF (ISTS.GT.0.AND.
     .        ISTS.LE.NLIM+NSTSI) THEN
C  ON NON DEFAULT SURFACE (ADD. OR STD.) ISTS=INMTI(IPOLGN,NRCELL)
        IF (ILIIN(ISTS) == 0) THEN
          NINCX=NCHBAR(IPOLGN,NRCELL)-NRCELL
        ELSE
          NINCX=SIGN(1,INMTI(IPOLGN,NRCELL))
        END IF
      ELSE
        GOTO 9999
      ENDIF
      IF (NLTRC) WRITE (iunout,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN A',
     .                         NRCELL,MRSURF,PT,NINCX,IPOLGN

      RETURN
      END IF

C PARTICLE OUTSIDE STANDARD MESH IN CELL NACELL, NRCELL=0 

      IF (.NOT.ALLOCATED(ITRINO)) THEN
        NSTS_CELL = COUNT(INMTI(1:3,1:NTRII) .NE. 0)
        ALLOCATE (ITRINO(NSTS_CELL))
        ALLOCATE (ISIDNO(NSTS_CELL))
        ISTS_CELL = 0
        DO ITRI=1,NTRII
          DO J=1,3
            IF (INMTI(J,ITRI) .NE. 0) THEN
              ISTS_CELL = ISTS_CELL + 1
              ITRINO(ISTS_CELL) = ITRI
              ISIDNO(ISTS_CELL) = J
            END IF
          END DO
        END DO
      END IF

      XX = X0
      YY = Y0
      PT=1.E30_DP
      ISTS_CELL = 0
      IF (NLSRFX) THEN
        NLSRFX=.FALSE.
        MPOL=MRSURF
        JPOL=IPOLG
      ELSE
        MPOL=0
        JPOL=0
      END IF
      MRSURF=0

      DO I=1, NSTS_CELL
        IZELL=ITRINO(I)
        J=ISIDNO(I)
        IF ((IZELL == MPOL) .AND. (J == JPOL)) CYCLE
        RICHTX = VTRIX(J,IZELL)
        RICHTY = VTRIY(J,IZELL)
        AX = XTRIAN(NECKE(J,IZELL))-XX
        AY = YTRIAN(NECKE(J,IZELL))-YY
        V  = (AX*VELY-AY*VELX)/(RICHTY*VELX-RICHTX*VELY+EPS60)
        IF (V.GE.0..AND.V.LE.1.) THEN
          IF (ABS(VELX).GT.ABS(VELY)) THEN
            T=(AX+V*RICHTX)/VELX
          ELSE
            T=(AY+V*RICHTY)/VELY
          ENDIF
          IF (NLTRC)
     .      WRITE (IUNOUT,*) ' INTERSECTION WITH TRIANGLE ',IZELL, 
     .                       ' SIDE ',J,' FOUND, PT,T = ',PT,T
          IF ((T .GT. 0.) .AND. ( T < PT)) THEN
            PT = T
            ISTS_CELL = I
          ENDIF
        ENDIF
      END DO

      IF (ISTS_CELL > 0) THEN
!  REENTRY FOUND
        NLSRFX=.TRUE.
        NJUMP=1
        MRSURF = ITRINO(ISTS_CELL)
        IPOLGN = ISIDNO(ISTS_CELL)
        NINCX = ITRINO(ISTS_CELL)
        IF (NLTRC) WRITE (iunout,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN B',
     .                                NRCELL,MRSURF,PT,NINCX,IPOLGN
      END IF

      RETURN
C
10000 CONTINUE

      IF (LEVGEO > 5) GOTO 12000
C
C     TETRAHEDRONS
C
      TIMTET = 1.E30
      NTIMT = 0
      IF (NRCELL .NE. 0) THEN
        IF (NJUMP.EQ.0) THEN
          IZELL = NRCELL
          IPOLGO=0
          IF (NLSRFX) THEN
            IPOLGO=IPOLG
          ENDIF
C       ELSEIF (NJUMP.EQ.1) THEN
        ENDIF
        IF (NLTRC) THEN
          WRITE (iunout,*) ' IZELL,IPOLGO ',IZELL,IPOLGO
          WRITE (iunout,*) ' X0,Y0,Z0 ',X0,Y0,Z0
          WRITE (iunout,*) ' VELX,VELY,VELZ ',VELX,VELY,VELZ
        ENDIF

        XX = X0
        YY = Y0
        ZZ = Z0
        TM=0.
        IF (IZELL.EQ.0) THEN
          WRITE (iunout,*) ' IZELL = 0 IN TIMER '
          CALL EXIT_OWN(1)

        END IF
C
C  CHECK TETRAHEDRON IZELL FOR INTERSECTION WITH TRAJECTORY
C
11020   CONTINUE
        DO J=1,4
          IF (J.NE.IPOLGO) THEN
            SIG1=SIGN(1,IRICH(1,J))
            I1=ABS(IRICH(1,J))
            SIG2=SIGN(1,IRICH(2,J))
            I2=ABS(IRICH(2,J))
c slmod begin - debug
            IF (izell.EQ.0) THEN 
              WRITE(6,*) 'NPANU =',npanu
              WRITE(6,*) 'IZELL =',izell
              WRITE(6,*) 'J     =',j
              WRITE(6,*) 'I1    =',i1
              WRITE(6,*) 'NRCELL=',nrcell
              WRITE(6,*) 'IPOLG =',ipolg
              WRITE(6,*) 'IPOLGO=',ipolgo
              PT=1.D30
              NINCX=0
              WRITE(iunout,*) 'ERROR (TIMER): Killing trajectory'
              RETURN
c              STOP 'HALTING THE CODE DUE TO TROUBLE'
            ENDIF
c slmod end
            PNORMI=RINCRC(J,IZELL)
            A(1:3,1) = (/ VTETX(I1,IZELL), VTETY(I1,IZELL),
     .                    VTETZ(I1,IZELL)/) * SIG1 * PNORMI
            A(1:3,2) = (/ VTETX(I2,IZELL), VTETY(I2,IZELL),
     .                    VTETZ(I2,IZELL)/) * SIG2 * PNORMI
            A(1:3,3) = (/ -VELX, -VELY, -VELZ /) * PNORMI
            B(1:3) = (/ XX-XTETRA(NTECK(ITSIDE(1,J),IZELL)),
     .                  YY-YTETRA(NTECK(ITSIDE(1,J),IZELL)),
     .                  ZZ-ZTETRA(NTECK(ITSIDE(1,J),IZELL)) /) * PNORMI
            DET = A(1,1) * A(2,2) * A(3,3)
     .          + A(1,2) * A(2,3) * A(3,1)
     .          + A(1,3) * A(2,1) * A(3,2)
     .          - A(3,1) * A(2,2) * A(1,3)
     .          - A(3,2) * A(2,3) * A(1,1)
     .          - A(3,3) * A(2,1) * A(1,2)
            IF (ABS(DET) < 1.D-5) CYCLE
            AB = A
            AB(:,1) = B
            XMU = ( AB(1,1) * AB(2,2) * AB(3,3)
     .            + AB(1,2) * AB(2,3) * AB(3,1)
     .            + AB(1,3) * AB(2,1) * AB(3,2)
     .            - AB(3,1) * AB(2,2) * AB(1,3)
     .            - AB(3,2) * AB(2,3) * AB(1,1)
     .            - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
            IF (XMU .GE.0.D0 .AND. XMU .LE.1.D0) THEN
              AB = A
              AB(:,2) = B
              XETA = ( AB(1,1) * AB(2,2) * AB(3,3)
     .               + AB(1,2) * AB(2,3) * AB(3,1)
     .               + AB(1,3) * AB(2,1) * AB(3,2)
     .               - AB(3,1) * AB(2,2) * AB(1,3)
     .               - AB(3,2) * AB(2,3) * AB(1,1)
     .               - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
              IF ((XETA.GE.0.D0 .AND. XETA.LE.1.D0) .AND.
     .            (XMU+XETA <= 1.D0)) THEN
                AB = A
                AB(:,3) = B
                T = ( AB(1,1) * AB(2,2) * AB(3,3)
     .              + AB(1,2) * AB(2,3) * AB(3,1)
     .              + AB(1,3) * AB(2,1) * AB(3,2)
     .              - AB(3,1) * AB(2,2) * AB(1,3)
     .              - AB(3,2) * AB(2,3) * AB(1,1)
     .              - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
                IF (T .GT. 0.) THEN
!            IF ((XMU .GE.0.D0 .AND. XMU .LE.1.D0) .AND.
!     .          (XETA.GE.0.D0 .AND. XETA.LE.1.D0) .AND.
!     .          (XMU+XETA <= 1.D0) .AND. (T .GT. 0.)) THEN
                  XX = XX + T * VELX
                  YY = YY + T * VELY
                  ZZ = ZZ + T * VELZ
                  if (nltrc) then
                     write(iunout,*) ' i1, i2 ',i1,i2
                     write(iunout,'(1x,a,3es14.7)') ' velxyz ',
     .                     velx,vely,velz
                     write(iunout,'(1x,a,3es14.7)') ' a11, a12, a13 ',
     .                     a(1,:)
                     write(iunout,'(1x,a,3es14.7)') ' a21, a22, a23 ',
     .                     a(2,:)
                     write(iunout,'(1x,a,3es14.7)') ' a31, a32, a33 ',
     .                     a(3,:)
                     write(iunout,'(1x,a,3es14.7)') ' b1, b2, b3 ',b(:)
                     write (iunout,'(1x,a,3es14.7)') ' det, xmu ',
     .                     det,xmu
                  end IF
                  TM = TM + T
                  TIMTET = TM
C  NUMMER DER GESCHNITTENEN SEITE DES ALTEN DREIECKS
                  NTIMT = J
                  IF (NLTRC) THEN
                    WRITE (iunout,*) 'IZELL,XX,YY,ZZ,J ',
     .                                IZELL,XX,YY,ZZ,J
                  ENDIF
C  ZELLENNUMMER DES NEUEN DREIECKS
                  IZELLO=IZELL
                  IPOLGOO=IPOLGO
                  IZELL = NTBAR(J,IZELLO)
                  PT = TM
                  IF (IZELL .EQ. 0) GOTO 11050
C  SEITENNUMMER DES NEUEN DREIECKS
                  IPOLGO = NTSEITE(J,IZELLO)
                  IF (NLTRC) THEN
                    WRITE(iunout,*) 'GEHE IN ZELLE ',IZELL,' TM= ',TM
                    WRITE(iunout,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
                  ENDIF
                  GOTO 11050
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        END DO
C
        IF (NLTRC) 
     .    WRITE (iunout,*) ' NO INTERSECTION FOUND IN TETRAHEDRON '
        PT=1.D30
        NINCX=0
C  DURCHSUCHE NACHBARDREIECK
        IF (NLSRFX) THEN
          IZELL=NTBAR(IPOLG,NRCELL)
          IPOLGO=NTSEITE(IPOLG,NRCELL)
c slmod begin
c          IF (IZELL.EQ.0) THEN
c            WRITE(6,*) 'IPOLG =',ipolg
c            WRITE(6,*) 'NRCELL=',nrcell
c            WRITE(6,*) 'IZELL =',izell
c            WRITE(6,*) 'IPOLGO=',ipolgo
c            WRITE (iunout,*) ' IZELL = 0 IN TIMER '
c            CALL EXIT_OWN(1)
c          END IF
c slmod end
          NRCELL=IZELL
          NLSRFX=.FALSE.
          GOTO 11020
        ENDIF
C  KEIN SCHNITTPUNKT GEFUNDEN UND NLSRFX=.FALSE.
        RETURN
C
C
C   THIS IS DONE FOR NJUMP=0 AND NJUMP=1
C
11050   CONTINUE
        NLSRFX=.FALSE.
        NJUMP=1
        MRSURF=NRCELL
        PT=TIMTET
        IPOLGN=NTIMT
        ISTS=ABS(INMTIT(IPOLGN,NRCELL))
        IF (ISTS.EQ.0) THEN
          NINCX=NTBAR(IPOLGN,NRCELL)-NRCELL
        ELSEIF (ISTS.GT.0.AND.
     .          ISTS.LE.NLIM+NSTSI) THEN
C  ON NON DEFAULT SURFACE (ADD. OR STD.) ISTS=INMTI(IPOLGN,NRCELL)
          IF (ILIIN(ISTS) == 0) THEN
            NINCX=NTBAR(IPOLGN,NRCELL)-NRCELL
          ELSE
            NINCX=SIGN(1,INMTIT(IPOLGN,NRCELL))
          END IF
        ELSE
          GOTO 9999
        ENDIF
        IF (NLTRC) WRITE (iunout,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN C',
     .                           NRCELL,MRSURF,PT,NINCX,IPOLGN

        RETURN

      ELSE
C
C  PARTICLE OUTSIDE STANDARD MESH, IN ADDITIONAL CELL NACELL
C  NRCELL=0
C
        IF (ITFRST == 0) THEN
          ITFRST = 1
          NTETSUR = COUNT(NTBAR(1:4,1:NTET) == 0)
          ALLOCATE (IDROB(NTETSUR,2))
          IOB = 0
          DO ITET = 1,NTET
            DO IS = 1,4
              IF (NTBAR(IS,ITET) == 0) THEN
                IOB = IOB + 1
                IDROB(IOB,:) = (/ ITET,IS /)
              END IF
            END DO
          END DO

          ALLOCATE (ISEEOB(0:NLIMI,NTETSUR/NBITS+1))
          ISEEOB = 0
          ISEEOB(0,:) = -1
          DO I=1,NTETSUR
            IT=IDROB(I,1)
            IS=IDROB(I,2)
            SIG1=SIGN(1,IRICH(1,IS))
            I1=ABS(IRICH(1,IS))
            SIG2=SIGN(1,IRICH(2,IS))
            I2=ABS(IRICH(2,IS))
            A1 = VTETX(I1,IT) * SIG1
            A2 = VTETY(I1,IT) * SIG1
            A3 = VTETZ(I1,IT) * SIG1
            B1 = VTETX(I2,IT) * SIG2
            B2 = VTETY(I2,IT) * SIG2
            B3 = VTETZ(I2,IT) * SIG2
            V1 = A2*B3 - A3*B2
            V2 = A3*B1 - A1*B3
            V3 = A1*B2 - A2*B1
            DO IL=1,NLIMI
              IF (RLB(IL) == 3.) THEN
                S1=(P1(1,IL)-XTETRA(NTECK(ITSIDE(1,IS),IT)))*V1 +
     .             (P1(2,IL)-YTETRA(NTECK(ITSIDE(1,IS),IT)))*V2 +
     .             (P1(3,IL)-ZTETRA(NTECK(ITSIDE(1,IS),IT)))*V3
                S2=(P2(1,IL)-XTETRA(NTECK(ITSIDE(1,IS),IT)))*V1 +
     .             (P2(2,IL)-YTETRA(NTECK(ITSIDE(1,IS),IT)))*V2 +
     .             (P2(3,IL)-ZTETRA(NTECK(ITSIDE(1,IS),IT)))*V3
                S3=(P3(1,IL)-XTETRA(NTECK(ITSIDE(1,IS),IT)))*V1 +
     .             (P3(2,IL)-YTETRA(NTECK(ITSIDE(1,IS),IT)))*V2 +
     .             (P3(3,IL)-ZTETRA(NTECK(ITSIDE(1,IS),IT)))*V3
                IF (ANY( (/ S1,S2,S3 /) > 0.D0)) THEN
                  CALL BITSET (ISEEOB,0,NLIMI,IL,I,1,NBITS)
                ELSE
                  IF (NOPTIM >= NTET) THEN
                    IF (NLIMPB >= NLIMPS) THEN
                      IGJUM3(ITET,IL) = 0
                    ELSE
                      CALL BITSET (IGJUM3,0,NOPTIM,ITET,IL,0,NBITS)
                    END IF
                  END IF
                END IF
              ELSE
                CALL BITSET (ISEEOB,0,NLIMI,IL,I,1,NBITS)
              END IF
            END DO
          END DO
        END IF

        PT=1.D30
        MRSURF=0

        TIMT=1.D30
        NTMZ=0
        NTMS=0
        DO IT=1,NTETSUR
          IF ((MRSURF == 0) .AND.
     .        .NOT.BITGET(ISEEOB,0,NLIMI,MSURF,IT,NBITS)) CYCLE
          I = IDROB(IT,1)
          J = IDROB(IT,2)
          IF (I == IZELLO)  CYCLE
          SIG1=SIGN(1,IRICH(1,J))
          I1=ABS(IRICH(1,J))
          SIG2=SIGN(1,IRICH(2,J))
          I2=ABS(IRICH(2,J))
          A(1:3,1) = (/ VTETX(I1,I), VTETY(I1,I),
     .                  VTETZ(I1,I)/) * SIG1
          A(1:3,2) = (/ VTETX(I2,I), VTETY(I2,I),
     .                  VTETZ(I2,I)/) * SIG2
          A(1:3,3) = (/ -VELX, -VELY, -VELZ /)
          B(1:3) = (/ X0-XTETRA(NTECK(ITSIDE(1,J),I)),
     .                Y0-YTETRA(NTECK(ITSIDE(1,J),I)),
     .                Z0-ZTETRA(NTECK(ITSIDE(1,J),I)) /)
          DET = SARRUS(A)
          IF (ABS(DET) < 1.D-5) CYCLE
          AB = A
          AB(:,1) = B
          XMU = SARRUS(AB)/DET
          AB = A
          AB(:,2) = B
          XETA = SARRUS(AB)/DET
          AB = A
          AB(:,3) = B
          T = SARRUS(AB)/DET
          IF ((XMU .GE.0.D0 .AND. XMU .LE.1.D0) .AND.
     .        (XETA.GE.0.D0 .AND. XETA.LE.1.D0) .AND.
     .        (XMU+XETA <= 1.D0) .AND. (T .GT. 0.)) THEN
            IF (T < TIMT) THEN
              TIMT = T
              NTMZ = I
              NTMS = J
            END IF
          END IF
        END DO
        PT = TIMT
        IZELLO = 0
        IPOLGOO = 0
        IF (NTMZ .NE. 0) THEN
C  SEITENNUMMER DES NEUEN DREIECKS
          IPOLGO = NTMS
          IF (NLTRC) THEN
            WRITE(iunout,*) 'OUTSIDE; GEHE IN ZELLE ',NTMZ,' TM= ',TM
            WRITE(iunout,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
          ENDIF
        END IF

        NLSRFX=.FALSE.
        NJUMP=1
        MRSURF=NTMZ
        IPOLGN=NTMS
        NINCX=NTMZ-NRCELL
        IF (NLTRC) WRITE (iunout,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN D',
     .                           NRCELL,MRSURF,PT,NINCX,IPOLGN
      END IF
      RETURN

12000 CONTINUE
C
C  GENERAL GEOMETRY OPTION: PROVIDE FLIGHT TIME IN CURRENT CELL
C
      IF (ICALL == 0) THEN
        ICALL = 1
        MXSF = MAXVAL(INUMP(1:NSTSI,1))
        IF (MXSF < NSURF) THEN
          NRMSRF = MXSF+1
        ELSE
          DO IRS=1,NSURF
            IF (ALL(INUMP(1:NSTSI,1).NE.IRS)) EXIT
          END DO
          NRMSRF = IRS
          IF (NRMSRF > NSURF) THEN
            WRITE (iunout,*) ' ERROR IN TIMER, NRMSRF WRONG '
            CALL EXIT_OWN(1)
          END IF
        END IF
      END IF
      CALL TIMUSR(NRCELL,X0,Y0,Z0,VELX,VELY,VELZ,NJUMP,
     .            NEWCEL,TIM,ICOS,IERR,NPANU,NLSRFX)
      IF (IERR.NE.0) GOTO 9999
      PT=TIM
      IF (NEWCEL.GT.0) THEN
        NINCX=NEWCEL-NRCELL
        MRSURF=NRMSRF
CPBDR   INMP1I(MRSURF,1,1)=0
      ELSEIF (NEWCEL.LT.0) THEN
        NINCX=ICOS
        MRSURF=INUMP(-NEWCEL,1)
CPBDR   INMP1I(MRSURF,1,1)=-NEWCEL
      ELSEIF (NEWCEL.EQ.0) THEN
        WRITE (iunout,*) 'NEWCEL=0, EXIT FROM MESH '
        CALL EXIT_OWN(1)
      ENDIF
      NLSRFX=.FALSE.
      RETURN

C
9999  CONTINUE
      WRITE (iunout,*) 'ERROR IN TIMER, EXIT CALLED AT 9999'
      WRITE (iunout,*) 'NPANU ', NPANU
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: timet.f
C
C
      SUBROUTINE TIMET (ZRAD)
C  INPUT
C
C   NRCELL:
C   NTCELL:
C   ZRAD  : DISTANCE (CM) TO THE NEXT RADIAL SURFACE OF 1D STANDARD MESH
C           OR TO NEXT ADDITIONAL SURFACE, TRAVELED IN RADIAL CELL NRCELL
C   PHI   :
C   X01   :
C   Z01   :
C
C  OUTPUT
C
C  IF NLTRZ AND NLTOR
C     NINCZ    :   DIRECTION IN Z-GRID
C     Z01      :
C     NTCELL   :
C     MTSURF   :
C  ELSEIF NLTRA
C    EITHER
C     ISRFCL<3  :   NO ROTATION , NNTCLL=0 , NO PARAMETERS CHANGED
C    OR
C     ISRFCL=3 :   ROTATION CLOCKWISE,
C                  STOP AND RESTART LATER AT MTSURF = NNTCLL
C     ISRFCL=3 :   ROTATION COUNTER CLOCKWISE
C                  STOP AND RESTART LATER AT MTSURF = NNTCLL+1
C
C     PHI      :
C     X01      :
C     Z01      :
C     ZRAD     :   REDUCED TO DISTANCE TO NEXT TOROIDAL SURFACE
C     NNTCLL   :   NEXT POSSIBLE CELL NUMBER IN TOROIDAL MESH
C                  IF PARTICLE TRAVELS REDUCED DISTANCE ZRAD
C     MTSURF   :   NEXT TOROIDAL SURFACE, IF ANY
C     NINCZ    :   DIRECTION IN Z-GRID
C     NINCX    :   RESET TO ZERO
C     MRSURF   :   RESET TO ZERO
C  ELSEIF NLTRP
C    TO BE WRITTEN
C
C  FOR THE TIME BEING: IF NLPOL, REDUCE PATH TO ONE TOROIDAL CELL
C
C
C     CALCULATE TIME SEGMENTS FOR 2D-PROFILES,
C     RADIALLY AND TOROIDALLY RESOLVED
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CGRID
      USE COMPRT
      USE COMSOU
      USE CLGIN

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: ZRAD
      REAL(DP) :: TTT, PHI0, X001, TO, AA, SUM, TU, XTO, EPSTST, XTU,
     .          F, FZ, ZRADS, BB, ZRD, DUM, Z001, X0TEST, Z0TEST, DZ,
     .          Y0TEST
      INTEGER :: ITT, ITEST, J2, NERR, IN, J, ICOU, MTTEST, ISTS,
     .           NZSAVE, INCZ, J1, IRSAVE, LEARCA
      INTEGER, SAVE :: MTSAVE=-1

      ZRADS=ZRAD
C
      IF (NLTRC) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'TIMET: ZRAD,NRCELL,NTCELL,MTSURF '
        WRITE (iunout,*) '      ',ZRAD,NRCELL,NTCELL,MTSURF
        WRITE (iunout,*) 'INITIAL: X0,X01,Z01,PHI ',
     .                             X0,X01,Z01,PHI/DEGRAD
      ENDIF
C
      IF (.NOT.NLTRZ) GOTO 1000
C
C  CYLINDRICAL OR CARTHESIAN CO-ORDINATE SYSTEM
C
C  PARTICLE OUTSIDE STANDARD MESH ?
C
      IF (NRCELL.EQ.0) THEN
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT TOROIDAL SURFACES
C
        NCOUT=0
        ZRD=ZRAD
        BB=VELZ+EPS60
        NZSAVE=1
        IF (VELZ.LT.0.D0) NZSAVE=-1
        DO 2903 ISTS=1,NSTSI
          MTTEST=INUMP(ISTS,3)
          IF (MTTEST.NE.0) THEN
C  TEST TOROIDAL SURFACE NO. MTTEST FOR REENTRY
C  TIME FROM Z01 TO ZSURF
            DZ=ZSURF(MTTEST)-Z01
            F=DZ/BB
            IF (NLTRC) WRITE (iunout,*) 'MTTEST,F,DZ ',MTTEST,F,DZ
            IF (F.LE.ZRD.AND.F.GT.0.D0) THEN
              X0TEST=X0+VELX*F
              Y0TEST=Y0+VELY*F
C  IS THIS RE-ENTRY INSIDE THE STANDARD GRID REGION?
              IF (LEVGEO.EQ.1) THEN
                IF (X0TEST.GE.RSURF(1).AND.X0TEST.LE.RSURF(NR1ST)) THEN
                  IRSAVE=LEARCA(X0TEST,RSURF,1,NR1ST,1,'TIMET 1    ')
                  NCOUT=1
                  JUPC(NCOUT)=1
                  MTSAVE=MTTEST
                  ZRD=F
                ENDIF
              ELSEIF (LEVGEO.GT.1) THEN
                WRITE (iunout,*) 'ERROR FROM TIMET: '
                WRITE (iunout,*) 
     .            'RE-ENTRY THROUGH TOROIDAL SURFACE NOT '
                WRITE (iunout,*) 'READY FOR LEVGEO.GT.1 '
                CALL EXIT_OWN(1)
              ENDIF
            ENDIF
          ENDIF
2903    CONTINUE
        IF (NCOUT.GT.0.AND.MTSAVE.NE.MTSURF) THEN
C  REENTRY FOUND, REDUCE ZRAD TO F
C  NCOUT=1 AT THIS POINT
          NCOUT=1
          MTSURF=MTSAVE
          NINCZ=NZSAVE
          IRCELL=IRSAVE
          ZRAD=ZRD
          ISRFCL=0
          BLPD(NCOUT)=ZRAD
          MRSURF=0
          MPSURF=0
          MASURF=0
          NINCX=0
          NINCY=0
          IF (NLTRC) THEN
            WRITE (iunout,*) 'REENTRY FOUND, MTSURF,ZRAD = ',MTSURF,ZRAD
            WRITE (iunout,*) 'IRCELL ',IRCELL
          ENDIF
          Z01=Z01+VELZ*ZRD
          NTCELL=KUPC(NCOUT)
          ITCELL=NTCELL
          GOTO 5000
        ELSE
C  NO REENTRY FOUND
          NCOUT=1
          KUPC(1)=1
          BLPD(1)=ZRAD
          IF (MRSURF.GT.0) THEN
C  CHECK VALID RANGE ON MRSURF
            Z0TEST=Z00+VELZ*ZRAD
            IF (NLTRC) WRITE (iunout,*) 
     .        'CHECK VALID RANGE: Z0TEST ',Z0TEST
            IF (Z0TEST.GE.ZSURF(1).AND.Z0TEST.LE.ZSURF(NT3RD)) THEN
              ITCELL=LEARCA(Z0TEST,ZSURF,1,NT3RD,1,'TIMET 2     ')
            ELSE
              MRSURF=0
              NINCX=0
            ENDIF
          ENDIF
          MTSURF=0
          NINCZ=0
          IF (NLTRC) THEN
            WRITE (iunout,*) 'NO REENTRY INTO TOROIDAL GRID FOUND '
          ENDIF
          Z01=Z01+ZRD*VELZ
          GOTO 5000
        ENDIF
C
      ENDIF
C
C  PARTICLE IN STANDARD MESH, RADIAL CELL NRCELL
C
2900  CONTINUE
      Z001=Z01+VELZ*ZRAD
C
      DUM=0.
      NCOUT=1
C
C  J1: CELL INDEX
C  J2  SURFACE INDEX
      J1=NTCELL
      IF (VELZ.LT.0.) THEN
        INCZ=0
        NINCZ=-1
        IF (NLSRFZ) J1=MTSURF-1
      ELSE
        INCZ=1
        NINCZ=1
        IF (NLSRFZ) J1=MTSURF
      ENDIF
      J2=J1+INCZ
C
      NLSRFZ=.FALSE.
C
10    CONTINUE
      IF (J2.LE.0.OR.J2.GT.NT3RD) GOTO 990
C  TIME FROM Z01 TO ZSURF
      DZ=(ZSURF(J2)-Z01)
!pb reduce zrad if the trajectory comes too close to the next standard surface
!pb without hitting it
      IF (ABS(ABS(DZ/ZRAD)-1._DP) <= EPS10) THEN
        ZRAD=(1._DP-EPS6)*ZRAD
        ZRADS=ZRAD
      END IF
      BB=VELZ+EPS60
      F=DZ/BB
      IF (F.LE.ZRAD) THEN
        KUPC(NCOUT)=J1
        BLPD(NCOUT)=F-DUM
        DUM=F
C  STOP HISTORY AT NON DEFAULT STANDARD SURFACE J2
        ITEST=INMP3I(0,0,J2)
        IN=ITEST+NLIM
!pb        IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
        IF ((ITEST.NE.0.AND.ILIIN(IN).NE.0).or.(ityp==3)) THEN
          ZRAD=F
          ISRFCL=0
          NTCELL=J1
          MTSURF=J2
          MRSURF=0
          MASURF=0
          NINCX=0
          IPOLGN=0
          ITCELL=NTCELL
          Z01=ZSURF(J2)
          IF (NLTRC) WRITE (iunout,*) 'STOP AT MTSURF ',MTSURF
          GOTO 5000
        ENDIF
C  NEXT CELL
        J1=J1+NINCZ
        J2=J1+INCZ
        NCOUT=NCOUT+1
        GOTO 10
      ENDIF
C
C  LAST CELL
C
100   CONTINUE
      NTCELL=J1
      MTSURF=0
      NINCZ=0
      KUPC(NCOUT)=J1
      BLPD(NCOUT)=ZRAD-DUM
      ITCELL=NTCELL
      Z01=Z001
C
      GOTO 5000
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN TIMET, Z SURFACE INDEX OUT OF RANGE  '
      WRITE (iunout,*) 'NPANU,Z0,Z01,ZRAD,VELZ,NTCELL '
      WRITE (iunout,*)  NPANU,Z0,Z01,ZRAD,VELZ,NTCELL
      WEIGHT=0.
      ZRAD=-1.
      GOTO 5000
C
C  CYLINDER APPROXIMATION FINISHED
C
1000  CONTINUE
C
      IF (.NOT.NLTRA) GOTO 5000
C
C  DISCRETE TOROIDAL APPROXIMATION, ONLY ONE STEP AT A TIME
C
      NERR=0
      NCOUT=1
      KUPC(1)=NTCELL
C
C     IF (NLSRFZ) THEN ....
      NLSRFZ=.FALSE.
C
1010  CONTINUE
C  PHI0 IS THE PHI AT THE CENTER OF THE CURRENT TOROIDAL CELL
      PHI0=PHI-ATAN2(Z01,X01)
C
      TTT=Z01/(X01*TANAL)
      IF (ABS(TTT).GT.1.+EPS10) THEN
        WRITE (iunout,*) 'NPANU ',NPANU
        WRITE (iunout,*) 'X01,Z01 OUT OF RANGE IN TIMET'
        WRITE (iunout,*) X01,Z01,TTT
        WRITE (iunout,*) 'TRY TO KILL PARTICLE ASAP '
        ZRAD=-1._DP
        RETURN
      ENDIF
C
      Z001=Z01+ZRAD*VELZ
      X001=X01+ZRAD*VELX
      IF (ZRAD.LT.1.D30.AND.X01*X001.GT.EPS10) THEN
        ITT=IDINT(REAL(Z001/(X001*TANAL),KIND(1.D0)))
        IF (NLTRC) WRITE (iunout,*) 'TIMET 1 ',X01,Z01,X001,Z001,ITT
      ELSE
        TO=(Z01-X01*TANAL)/(TANAL*VELX-VELZ)
        XTO=X01+TO*VELX
        TU=(Z01+X01*TANAL)/(-TANAL*VELX-VELZ)
        XTU=X01+TU*VELX
        IF (NLTRC) WRITE (iunout,*) 'TU,TO ',TU,TO,XTU,XTO
        EPSTST=EPS10*TANAL
        IF (XTO.GT.0..AND.TO.GT.EPSTST) THEN
          ITT=1
        ELSEIF (XTU.GT.0..AND.TU.GT.EPSTST) THEN
          ITT=-1
        ELSE
          ITT=0
          Z001=Z01
          X001=X01
        ENDIF
        IF (NLTRC) WRITE (iunout,*) 'TIMET 2 ',X01,Z01,ITT
      ENDIF
C
      IF (ITT) 1100,1200,1300
C
C  NO INTERSECTION WITH TOROIDAL SURFACE
C
1200  CONTINUE
      MTSURF=0
      NNTCLL=IPERID
      Z01=Z001
      X01=X001
C  IN CASE ZRAD=1.D30, IS NEXT STATEMENT IS NONSENSE, BUT CORRECTED
C                      FOR IN SUBR. STDCOL, WHICH MUST BE CALLED NEXT
C                      FOR A POLOIDAL SURFACE (OTHERWISE: ERROR EXIT)
      PHI=PHI0+ATAN2(Z01,X01)
      BLPD(1)=ZRAD
      IF (NLTRC) WRITE (iunout,*) 'FINAL 1: X01,Z01,PHI ',
     .                                      X01,Z01,PHI/DEGRAD
      GOTO 5000
C
C  PARTICLE LEAVES CELL IN POSITIVE DIRECTION, REDUCE ZRAD
C
1300  CONTINUE
C
C  TIME TO REACH CELL SURFACE: F
C
      AA=Z01-X01*TANAL
      BB=(TANAL*VELX-VELZ)+EPS60
      F=AA/BB
      IF (F.LE.0.) THEN
C  PARTICLE ACCIDENTALLY ON A TOROIDAL SURFACE?
        IF (ABS(AA).LE.EPS10.AND.NERR.LE.1) THEN
          IF (NLTRC) WRITE (iunout,*) 'TRY AGAIN IN TIMET'
          Z01=Z01-VELZ*EPS10
          X01=X01-VELX*EPS10
          NERR=NERR+1
          GOTO 1010
        ENDIF
        ZRAD=1.D30
        GOTO 9998
      ENDIF
      X01=X01+F*VELX
      Z01=TANAL*X01
      PHI=PHI0+ZHALF
      ZRAD=F
C
      ISRFCL=3
      BLPD(1)=ZRAD
C
      NINCX=0
      NINCZ=1
      IPOLGN=0
      MRSURF=0
      MASURF=0
      NNTCLL=NTCELL+1
      IF (.NOT.NLTOR) NNTCLL=IPERID+1
      MTSURF=NNTCLL
C  ENFORCED PERIODICITY, UNLESS NON-DEFAULT STANDARD SURFACE
      ISTS=0
      IF (NLTOR) ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
      IF (ISTS.EQ.0) THEN
        IF (NNTCLL.GE.NTTRA) NNTCLL=1
        MTSURF=NNTCLL
      ELSE
C  NO AUTOMATIC PERIODICITY IN SUBR. TORCOL. CALL STDCOL FOR NON DEF. SURF.
        ISRFCL=0
      ENDIF
C
      IF (NLTRC) THEN
        WRITE (iunout,*) 'FINAL 2: X01,Z01,PHI ',X01,Z01,PHI/DEGRAD
        WRITE (iunout,*) 'ZRAD,ISTS,ISRFCL,IPERID ',
     .                    ZRAD,ISTS,ISRFCL,IPERID
      ENDIF
      GOTO 5000
C
C  PARTICLE LEAVES CELL IN NEGATIVE DIRECTION, REDUCE ZRAD
C
1100  CONTINUE
C
C  TIME TO REACH CELL SURFACE
C
      AA=Z01+X01*TANAL
      BB=-TANAL*VELX-VELZ+EPS60
      F=AA/BB
      IF (F.LE.0.) THEN
C  PARTICLE ACCIDENTALLY ON A TOROIDAL SURFACE?
        IF (ABS(AA).LE.EPS10.AND.NERR.LE.1) THEN
          IF (NLTRC) WRITE (iunout,*) 'TRY AGAIN IN TIMET'
          Z01=Z01-VELZ*EPS10
          X01=X01-VELX*EPS10
          NERR=NERR+1
          GOTO 1010
        ENDIF
        ZRAD=1.D30
        GOTO 9998
      ENDIF
      X01=X01+F*VELX
      Z01=-TANAL*X01
      PHI=PHI0-ZHALF
      ZRAD=F
C
      ISRFCL=3
      BLPD(1)=ZRAD
C
      NINCX=0
      NINCZ=-1
      IPOLGN=0
      MRSURF=0
      MASURF=0
      NNTCLL=NTCELL-1
      IF (.NOT.NLTOR) NNTCLL=IPERID-1
      MTSURF=NNTCLL+1
C  ENFORCED PERIODICITY, UNLESS NON-DEFAULT STANDARD SURFACE
      ISTS=0
      IF (NLTOR) ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
      IF (ISTS.EQ.0) THEN
        IF (NNTCLL.LE.0) NNTCLL=NTTRAM
        MTSURF=NNTCLL+1
      ELSE
C  NO AUTOMATIC PERIODICITY IN SUBR. TORCOL. CALL STDCOL FOR NON DEF. SURF.
        ISRFCL=0
      ENDIF
C
      IF (NLTRC) THEN
        WRITE (iunout,*) 'FINAL 3: X01,Z01,PHI ',X01,Z01,PHI/DEGRAD
        WRITE (iunout,*) 'ZRAD,ISTS,ISRFCL,IPERID ',
     .                    ZRAD,ISTS,ISRFCL,IPERID
      ENDIF
      GOTO 5000
C
C  DISCRETE TOROIDAL APPROXIMATION FINISHED
C
C
5000  CONTINUE
      IF (NLTRC) WRITE (iunout,*) 'NCOUT= ',NCOUT
      DO 5100 J=1,NCOUT
        CLPD(J)=BLPD(J)
        NUPC(J)=(KUPC(J)-1)*NP2T3
        NCOUNT(J)=KUPC(J)
        IF (CLPD(J).LE.0..OR.KUPC(J).LE.0.OR.
     .      (KUPC(J).GE.NT3RD.AND.NLTOR)) THEN
          WRITE (iunout,*) 'ERROR DETECTED IN TIMET '
          WRITE (iunout,*) 'NPANU,J,BLPD,KUPC ',NPANU,J,BLPD(J),KUPC(J)
        ENDIF
        IF (NLTRC) THEN
          WRITE (iunout,*) 'TIMET: J,BLPD,NUPC,NCOUNT ',
     .                        J,BLPD(J),NUPC(J),NCOUNT(J)
        ENDIF
5100  CONTINUE
      IF (NLTRC) THEN
        WRITE (iunout,*) 'MTSURF,NLSRFZ,NINCZ,IRCELL,IPCELL '
        WRITE (iunout,*)  MTSURF,NLSRFZ,NINCZ,IRCELL,IPCELL
        IF (NLTOR) 
     .    WRITE (iunout,*) 'INMP3I ',INMP3I(IRCELL,IPCELL,MTSURF)
      ENDIF
C
      NCOU=NCOUT
C
      SUM=0.
      DO 5110 ICOU=1,NCOU
        SUM=SUM+CLPD(ICOU)
C       WRITE (iunout,*) 'ICOU,NCOUNT,CLPD ',
C    .                    ICOU,NCOUNT(ICOU),CLPD(ICOU)
5110  CONTINUE
      IF (MTSURF.EQ.0.AND.ABS(SUM-ZRADS).GT.EPS10) THEN
        WRITE (iunout,*) 'ERROR IN TIMET: NPANU,SUM,ZRADS ',
     .                               NPANU,SUM,ZRADS
        WRITE (iunout,*) 'TRY TO KILL PARTICLE ASAP '
        SUM=-1.0D0
      ENDIF
      ZRAD=SUM
C
      RETURN
C
9998  CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN TIMET. RETURN ZRAD=1.D30'
      WRITE (iunout,*) 'NPANU,AA,BB ',NPANU,AA,BB
      RETURN
9999  CONTINUE
      WRITE (iunout,*) 'INVALID OPTION IN TIMET. EXIT CALLED '
      CALL EXIT_OWN(1)
      END
