C
C
      SUBROUTINE TIMEA
C
C   1 ST INTERSECTION OF THE RAY X+T*VX,Y+T*VY,Z+T*VZ WITH ONE OF THE
C   ADDITIONAL SURFACES, DEFINED BY 2.ND ORDER EQUATIONS
C   IT IS ALSO CHECKED, WETHER THIS INTERSECTION TAKES PLACE INSIDE THE
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
      LOGICAL :: LMTSRF(NLIMPS)
      SAVE
C
      ENTRY TIMEA0
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
          A0LM(J)=DZ3*P1(2,J)-DY2*P1(3,J)
          A1LM(J)=0.
          A2LM(J)=-DZ3
          A3LM(J)=DY2
          YLIMS1(1,J)=MIN(P1(2,J),P2(2,J))
          YLIMS2(1,J)=MAX(P1(2,J),P2(2,J))
          ZLIMS1(1,J)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(1,J)=MAX(P1(3,J),P2(3,J))
!          
!          write (6,'(a,6es12.4)') 'xylims1,2 ',
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
          A0LM(J)=DZ3*P1(1,J)-DX1*P1(3,J)
          A1LM(J)=-DZ3
          A2LM(J)=0.
          A3LM(J)=DX1
          XLIMS1(1,J)=MIN(P1(1,J),P2(1,J))
          XLIMS2(1,J)=MAX(P1(1,J),P2(1,J))
          ZLIMS1(1,J)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(1,J)=MAX(P1(3,J),P2(3,J))
!          
!          write (6,'(a,6es12.4)') 'xylims1,2 ',
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
          A0LM(J)=DY2*P1(1,J)-DX1*P1(2,J)
          A1LM(J)=-DY2
          A2LM(J)=DX1
          A3LM(J)=0.
          XLIMS1(1,J)=MIN(P1(1,J),P2(1,J))
          XLIMS2(1,J)=MAX(P1(1,J),P2(1,J))
          YLIMS1(1,J)=MIN(P1(2,J),P2(2,J))
          YLIMS2(1,J)=MAX(P1(2,J),P2(2,J))
!          
!          write (6,'(a,6es12.4)') 'xylims1,2 ',
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
          WRITE (6,*) 'WARNING FROM TIMEA0, TEST FOR 4.TH POINT:'
          WRITE (6,*) 'SF. NUMBER, TEST= ',J,TEST4
        ELSEIF (ABS(TEST5).GT.EPS10) THEN
          WRITE (6,*) 'WARNING FROM TIMEA0, TEST FOR 5.TH POINT:'
          WRITE (6,*) 'SF. NUMBER, TEST= ',J,TEST5
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
          WRITE (6,*) 'EXIT FROM TIMEA0: SURFACE NO J IS OPERATING'
          WRITE (6,*) 'A SWITCH BUT IS NOT SEEN BY HISTORY'
          WRITE (6,*) 'J= ',J
          CALL EXIT_OWN(1)
        ENDIF
C
        IF (ILSWCH(J).NE.0.AND.ILIIN(J).GT.0) THEN
          DO ISPZ=1,NSPTOT
            IF (TRANSP(ISPZ,1,J).NE.0.D0.OR.TRANSP(ISPZ,2,J).NE.0.D0)
     .      GOTO 95
          ENDDO
          GOTO 96
95        WRITE (6,*) 'EXIT FROM TIMEA0: SURFACE NO J IS OPERATING'
          WRITE (6,*) 'A SWITCH BUT IS SOMETIMES TRANSPARENT AND '
          WRITE (6,*) 'SOMETIMES REFLECTING (SEMI-TRANSPARENCY OPTION)'
          WRITE (6,*) 'POSSIBLE FIXES: MANUAL, CHAPTER 2, SECTION 6 '
          WRITE (6,*) 'J= ',J
          CALL EXIT_OWN(1)
        ENDIF
C
96      GOTO 94
C
98      WRITE (6,*) 'WARNING FROM TIMEA0, SURFACE NO J IS TURNED OF'
        WRITE (6,*) 'BECAUSE ILL DEFINED '
        WRITE (6,*) 'J= ',J
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
C     IF (NLTRC) THEN
C       CALL LEER(1)
C       WRITE (6,*) 'TIMEA ,X,Y,Z,T ',X,Y,Z,T
C       WRITE (6,*) '      VX,VY,VZ,V ',VX,VY,VZ,V
C       IF (NLTRA) WRITE (6,*) 'MSURF,NNTCL ',MSURF,NNTCL
C       IF (.NOT.NLTRA) WRITE (6,*) 'MSURF ',MSURF
C     ENDIF
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
C           WRITE (6,*) 'J,ILTOR(J),X,Z,VX,VZ ',J,ILTOR(J),X,Z,VX,VZ
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
C           WRITE (6,*) 'J,ILTOR(J),X,Z,VX,VZ ',J,ILTOR(J),X,Z,VX,VZ
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
C       IF (NLTRC) WRITE (6,*) 'TIMEA, J,TMI,TMA ',J,TMI,TMA
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
C       IF (NLTRC) WRITE (6,*) 'TIMEA AT 40, J,TMA,A1,A2 ',J,TMA,A1,A2
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
C       IF (NLTRC) WRITE (6,*) 'TIMEA AT 50, J,TMA,A2,A3 ',J,TMA,A2,A3
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
C       IF (NLTRC) THEN
C         WRITE (6,*) 'TIMEA, TL,XR,YR,ZR,NNR,MASURF,SG '
C         WRITE (6,*)         TL,XR,YR,ZR,NNR,MASURF,SG
C       ENDIF
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
C     IF (NLTRC) WRITE (6,*) 'NOT RETURNED FROM TIMEA, OTHER LOOP '
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
