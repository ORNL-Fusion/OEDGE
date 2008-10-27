      SUBROUTINE SETFIT(TRCSUR)
c new option, 2002: ilfit.lt.0   not fully tested
c                   connecting surface "ie" not necessarily rlb=1, or 1.5
c                   anymore.
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE CCONA
      USE CLGIN

      IMPLICIT NONE
C
      LOGICAL, INTENT(IN) :: TRCSUR
      REAL(DP) :: XS(8), YS(8), D, DST(8,2), XA(14), XB(14)
      REAL(DP) :: XT, SRAD, RAD, YMIN, YMAX, B, AH, YT, 
     .          DX1, DY2, DZ3, XNORM, AT, A3, A4, A2, XR, XL, A5,
     .          A1, A0, YR, YL, B0, B1, B2, B3, B4, XMAX, XMIN, XP1,
     .          XP2, XP3, XP4, YP1, YP2, YP3, YP4, B5
      INTEGER :: IEQ(2), ILFT, ILFTS
      INTEGER :: IMIN1, IMIN2, K, IS, JUM, IPNT1, IPNT2, I, IE, J
      LOGICAL :: LINFX,LINFY,LINFZ, TWOPOINT
C
      DO 1 I=1,NLIMI
C  IS ILFIT OPTION IN USE?
        IF (ILFIT(I).EQ.0.OR.IGJUM0(I).NE.0) GOTO 1
C  SELECT THE SURFACE NUMBERS OF THE NEIGHBORING SURFACES
        IF (ILFIT(I).LT.0) THEN
          ILFT=-ILFIT(I)
          ILFTS=-1
        ELSE
          ILFT=ILFIT(I)
          ILFTS=1
        ENDIF
        IEQ(1)=ILFT/1000
        IEQ(2)=ILFT-IEQ(1)*1000
C  SURFACE I MUST BE GIVEN BY TWO POINT OPTION WITH ONE
C  IGNORABLE CO-ORDINATE
C  USE THIRD POINT FOR IDENTIFICATION OF 2-POINT INPUT OPTION
C  BECAUSE RLB HAS ALREADY BEEN OVERWRITTEN
        IF (P3(1,I).LT.1.D50.AND.P3(2,I).LT.1.D50.AND.
     .      P3(3,I).LT.1.D50) GOTO 991
C  WHICH CO-ORDINATE OF SURFACE NO. I IS IGNORABLE?
        LINFX=P3(1,I).GT.1.D50
        LINFY=P3(2,I).GT.1.D50
        LINFZ=P3(3,I).GT.1.D50
        IPNT1=0
        IPNT2=0
        DO 2 J=1,2
          IE=IEQ(J)
          IF (IE.EQ.0) GOTO 2
          IF (IGJUM0(IE).NE.0) GOTO 991
C  CONNECT SURFACE I WITH SURFACE NO. IE
C  SURFACE IE MUST BE GIVEN WITH SAME IGNORABLE CO-ORDINATE AS SURFACE I
C
C  IDENTIFY IGNORABLE CO-ORDINATE OF SURFACE IE
          IF (LINFX) THEN
            IF (A1LM(IE).NE.0..OR.A4LM(IE).NE.0..OR.
     .          A7LM(IE).NE.0..OR.A8LM(IE).NE.0.D0) GOTO 5
C  X IS IGNORABLE, IN BOTH SURFACES: I AND IE
            A0=A0LM(IE)
            A1=A2LM(IE)
            A2=A3LM(IE)
            A3=A5LM(IE)
            A4=A6LM(IE)
            A5=A9LM(IE)
            XL=YLIMS1(1,IE)
            XR=YLIMS2(1,IE)
            YL=ZLIMS1(1,IE)
            YR=ZLIMS2(1,IE)

            XP1=P1(2,I)
            XP2=P2(2,I)
            YP1=P1(3,I)
            YP2=P2(3,I)
C
            XP3=P1(2,IE)
            XP4=P2(2,IE)
            YP3=P1(3,IE)
            YP4=P2(3,IE)
C                      
            GOTO 100
          ENDIF
5         CONTINUE
          LINFX=.FALSE.
          IF (LINFY) THEN
            IF (A2LM(IE).NE.0..OR.A5LM(IE).NE.0..OR.
     .          A7LM(IE).NE.0..OR.A9LM(IE).NE.0.D0) GOTO 6
C  Y IS IGNORABLE, IN BOTH SURFACES: I AND IE
            A0=A0LM(IE)
            A1=A1LM(IE)
            A2=A3LM(IE)
            A3=A4LM(IE)
            A4=A6LM(IE)
            A5=A8LM(IE)
            XL=XLIMS1(1,IE)
            XR=XLIMS2(1,IE)
            YL=ZLIMS1(1,IE)
            YR=ZLIMS2(1,IE)

            XP1=P1(1,I)
            XP2=P2(1,I)
            YP1=P1(3,I)
            YP2=P2(3,I)
C     
            XP3=P1(1,IE)
            XP4=P2(1,IE)
            YP3=P1(3,IE)
            YP4=P2(3,IE)
C
            GOTO 100
          ENDIF
6         CONTINUE
          LINFY=.FALSE.
          IF (LINFZ) THEN
            IF (A3LM(IE).NE.0..OR.A6LM(IE).NE.0..OR.
     .          A8LM(IE).NE.0..OR.A9LM(IE).NE.0.D0) GOTO 7
C  Z IS IGNORABLE, IN BOTH SURFACES: I AND IE
            A0=A0LM(IE)
            A1=A1LM(IE)
            A2=A2LM(IE)
            A3=A4LM(IE)
            A4=A5LM(IE)
            A5=A7LM(IE)
            XL=XLIMS1(1,IE)
            XR=XLIMS2(1,IE)
            YL=YLIMS1(1,IE)
            YR=YLIMS2(1,IE)

            XP1=P1(1,I)
            XP2=P2(1,I)
            YP1=P1(2,I)
            YP2=P2(2,I)
C
            XP3=P1(1,IE)
            XP4=P2(1,IE)
            YP3=P1(2,IE)
            YP4=P2(2,IE)
C
            GOTO 100
          ENDIF
7         CONTINUE
          LINFZ=.FALSE.
          GOTO 990
100       CONTINUE
C
C
C  IS SURFACE IE GIVEN BY TWO-POINT OPTION ?
        TWOPOINT=.FALSE.
        IF (P3(1,IE).GE.1.D50.OR.P3(2,IE).GE.1.D50.OR.
     .      P3(3,IE).GE.1.D50) TWOPOINT=.TRUE.

          IF ((RLB(IE).EQ.1..OR.RLB(IE).EQ.1.5).AND..NOT.TWOPOINT) THEN
C  TRY TO CONNECT SURFACE I TO POINTS ON BOUNDARY BOX OF SURFACE IE
          IS=0
C  SCHNITTPUNKTE VON IE MIT X=XL AND X=XR
          IF (ABS(A4).LE.EPS12) THEN
            YT=-(A0+A1*XL+A3*XL*XL)/(A2+A5*XL)
            IF (YT.GE.YL.AND.YT.LE.YR) THEN
              IS=IS+1
              XS(IS)=XL
              YS(IS)=YT
            ENDIF
            YT=-(A0+A1*XR+A3*XR*XR)/(A2+A5*XR)
            IF (YT.GE.YL.AND.YT.LE.YR) THEN
              IS=IS+1
              XS(IS)=XR
              YS(IS)=YT
            ENDIF
          ELSE
            AH=0.5*(A2+A5*XL)/A4
            B=(A0+A1*XL+A3*XL*XL)/A4
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              YT=-AH+SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XL
                YS(IS)=YT
              ENDIF
              YT=-AH-SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XL
                YS(IS)=YT
              ENDIF
            ENDIF
            AH=0.5*(A2+A5*XR)/A4
            B=(A0+A1*XR+A3*XR*XR)/A4
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              YT=-AH+SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XR
                YS(IS)=YT
              ENDIF
              YT=-AH-SRAD
              IF (YT.GE.YL.AND.YT.LE.YR) THEN
                IS=IS+1
                XS(IS)=XR
                YS(IS)=YT
              ENDIF
            ENDIF
          ENDIF
C  SCHNITTPUNKTE MIT Y=YL AND Y=YR
          IF (ABS(A3).LE.EPS12) THEN
            XT=-(A0+A2*YL+A4*YL*YL)/(A1+A5*YL)
            IF (XT.GE.XL.AND.XT.LE.XR) THEN
              IS=IS+1
              YS(IS)=YL
              XS(IS)=XT
            ENDIF
            XT=-(A0+A2*YR+A4*YR*YR)/(A1+A5*YR)
            IF (XT.GE.XL.AND.XT.LE.XR) THEN
              IS=IS+1
              YS(IS)=YR
              XS(IS)=XT
            ENDIF
          ELSE
            AH=0.5*(A1+A5*YL)/A3
            B=(A0+A2*YL+A4*YL*YL)/A3
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              XT=-AH+SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YL
                XS(IS)=XT
              ENDIF
              XT=-AH-SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YL
                XS(IS)=XT
              ENDIF
            ENDIF
            AH=0.5*(A1+A5*YR)/A3
            B=(A0+A2*YR+A4*YR*YR)/A3
            RAD=AH*AH-B
            IF (RAD.GE.0.D0) THEN
              SRAD=SQRT(RAD)
              XT=-AH+SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YR
                XS(IS)=XT
              ENDIF
              XT=-AH-SRAD
              IF (XT.GE.XL.AND.XT.LE.XR) THEN
                IS=IS+1
                YS(IS)=YR
                XS(IS)=XT
              ENDIF
            ENDIF
          ENDIF

          ELSEIF ((RLB(IE).EQ.1..OR.RLB(IE).EQ.1.5).AND.TWOPOINT) THEN
C  TRY TO CONNECT SURFACE I TO SURFACE IE
            IS=1
            CALL SCHNITP(XP1,YP1,XP2,YP2,XP3,YP3,XP4,YP4,XS(1),YS(1))
          ELSE
C  ALL OTHER RLB(IE) OPTIONS
            GOTO 991
          ENDIF
C
C  SELECT THE DISTANCES TO THE INTERSECTION POINTS
          DO 3 K=1,IS
            DST(K,1)=1.D60
            DST(K,2)=1.D60
            IF (IPNT1.EQ.0)
     .      DST(K,1)=SQRT((XS(K)-XP1)**2+(YS(K)-YP1)**2)
            IF (IPNT2.EQ.0)
     .      DST(K,2)=SQRT((XS(K)-XP2)**2+(YS(K)-YP2)**2)
3         CONTINUE
C
          IMIN1=1
          IMIN2=1
          DO 4 K=2,IS
            IF (DST(K,1).LT.DST(IMIN1,1)) IMIN1=K
            IF (DST(K,2).LT.DST(IMIN2,2)) IMIN2=K
4         CONTINUE
C
          IF (TRCSUR) THEN
            WRITE (6,*) 'SETFIT, I,IE ',I,IE
            DO 4711 K=1,IS
              WRITE (6,*) 'K,XS,YS,DIST1,DIST2 ',
     .                     K,XS(K),YS(K),DST(K,1),DST(K,2)
4711        CONTINUE
          ENDIF
C
C  SET THE SELECTED POINT
          IF (ILFTS.LT.0) THEN
C  EXCHANGE POINTS TO BE RESET: NOT THE CLOSEST, BUT THE OTHER ONE
            D=DST(IMIN1,1)
            DST(IMIN1,1)=DST(IMIN2,2)
            DST(IMIN2,2)=D
          ENDIF
          IF (DST(IMIN1,1).LT.DST(IMIN2,2)) THEN
            XP1=XS(IMIN1)
            YP1=YS(IMIN1)
            IPNT1=1
            IF (TRCSUR)
     .      WRITE (6,*) 'PNT.1 REPLACED BY ',XP1,YP1,' FOR SURF. ',I
          ELSE
            XP2=XS(IMIN2)
            YP2=YS(IMIN2)
            IPNT2=1
            IF (TRCSUR)
     .      WRITE (6,*) 'PNT.2 REPLACED BY ',XP2,YP2,' FOR SURF. ',I
          ENDIF
C
C  RESET THE POINTS ON THE ARRAYS
          IF (TRCSUR) WRITE (6,*) ' LINFXYZ ',LINFX,LINFY,LINFZ
          IF (LINFX) THEN
            P1(2,I)=XP1
            P1(3,I)=YP1
            P2(2,I)=XP2
            P2(3,I)=YP2
            DZ3=(P1(3,I)-P2(3,I))
            DY2=(P1(2,I)-P2(2,I))
            A0LM(I)=DZ3*P1(2,I)-DY2*P1(3,I)
            A1LM(I)=0.
            A2LM(I)=-DZ3
            A3LM(I)=DY2
            YLIMS1(1,I)=MIN(P1(2,I),P2(2,I))
            YLIMS2(1,I)=MAX(P1(2,I),P2(2,I))
            ZLIMS1(1,I)=MIN(P1(3,I),P2(3,I))
            ZLIMS2(1,I)=MAX(P1(3,I),P2(3,I))
            IF (P1(2,I).EQ.P2(2,I)) THEN
              YLIMS1(1,I)=YLIMS1(1,I)-0.1
              YLIMS2(1,I)=YLIMS2(1,I)+0.1
            ENDIF
            IF (P1(3,I).EQ.P2(3,I)) THEN
              ZLIMS1(1,I)=ZLIMS1(1,I)-0.1
              ZLIMS2(1,I)=ZLIMS2(1,I)+0.1
            ENDIF
          ELSEIF (LINFY) THEN
            P1(1,I)=XP1
            P1(3,I)=YP1
            P2(1,I)=XP2
            P2(3,I)=YP2
            DZ3=(P1(3,I)-P2(3,I))
            DX1=(P1(1,I)-P2(1,I))
            A0LM(I)=DZ3*P1(1,I)-DX1*P1(3,I)
            A1LM(I)=-DZ3
            A2LM(I)=0.
            A3LM(I)=DX1
            XLIMS1(1,I)=MIN(P1(1,I),P2(1,I))
            XLIMS2(1,I)=MAX(P1(1,I),P2(1,I))
            ZLIMS1(1,I)=MIN(P1(3,I),P2(3,I))
            ZLIMS2(1,I)=MAX(P1(3,I),P2(3,I))
            IF (P1(1,I).EQ.P2(1,I)) THEN
              XLIMS1(1,I)=XLIMS1(1,I)-0.1
              XLIMS2(1,I)=XLIMS2(1,I)+0.1
            ENDIF
            IF (P1(3,I).EQ.P2(3,I)) THEN
              ZLIMS1(1,I)=ZLIMS1(1,I)-0.1
              ZLIMS2(1,I)=ZLIMS2(1,I)+0.1
            ENDIF
          ELSEIF (LINFZ) THEN
            P1(1,I)=XP1
            P1(2,I)=YP1
            P2(1,I)=XP2
            P2(2,I)=YP2
            DY2=(P1(2,I)-P2(2,I))
            DX1=(P1(1,I)-P2(1,I))
            A0LM(I)=DY2*P1(1,I)-DX1*P1(2,I)
            A1LM(I)=-DY2
            A2LM(I)=DX1
            A3LM(I)=0.
            XLIMS1(1,I)=MIN(P1(1,I),P2(1,I))
            XLIMS2(1,I)=MAX(P1(1,I),P2(1,I))
            YLIMS1(1,I)=MIN(P1(2,I),P2(2,I))
            YLIMS2(1,I)=MAX(P1(2,I),P2(2,I))
            IF (P1(1,I).EQ.P2(1,I)) THEN
              XLIMS1(1,I)=XLIMS1(1,I)-0.1
              XLIMS2(1,I)=XLIMS2(1,I)+0.1
            ENDIF
            IF (P1(2,I).EQ.P2(2,I)) THEN
              YLIMS1(1,I)=YLIMS1(1,I)-0.1
              YLIMS2(1,I)=YLIMS2(1,I)+0.1
            ENDIF
          ENDIF
C
          AT=MAX(ABS(A1LM(I)),ABS(A2LM(I)),ABS(A3LM(I)))
          IF (ABS(A1LM(I)).EQ.AT) JUMLIM(I)=1
          IF (ABS(A2LM(I)).EQ.AT) JUMLIM(I)=2
          IF (ABS(A3LM(I)).EQ.AT) JUMLIM(I)=3
          XNORM=SQRT(A1LM(I)*A1LM(I)+A2LM(I)*A2LM(I)+A3LM(I)*A3LM(I))
          IF (XNORM.LE.EPS60) GOTO 993
          A0LM(I)=A0LM(I)/XNORM
          A1LM(I)=A1LM(I)/XNORM
          A2LM(I)=A2LM(I)/XNORM
          A3LM(I)=A3LM(I)/XNORM
          JUM=JUMLIM(I)
          GOTO (91,92,93),JUM
91          ALM(I)=-A0LM(I)/A1LM(I)
            BLM(I)=-A2LM(I)/A1LM(I)
            CLM(I)=-A3LM(I)/A1LM(I)
          GOTO 97
92          ALM(I)=-A0LM(I)/A2LM(I)
            BLM(I)=-A1LM(I)/A2LM(I)
            CLM(I)=-A3LM(I)/A2LM(I)
          GOTO 97
93          ALM(I)=-A0LM(I)/A3LM(I)
            BLM(I)=-A1LM(I)/A3LM(I)
            CLM(I)=-A2LM(I)/A3LM(I)
97        CONTINUE
C
          IF (TRCSUR) THEN
            WRITE (6,*) ' A0-A3 ',A0LM(I),A1LM(I),A2LM(I),A3LM(I)
            WRITE (6,*) ' XLIMS ',XLIMS1(1,I),XLIMS2(1,I)
            WRITE (6,*) ' YLIMS ',YLIMS1(1,I),YLIMS2(1,I)
            WRITE (6,*) ' ZLIMS ',ZLIMS1(1,I),ZLIMS2(1,I)
          ENDIF
C
C  TWO POINT OPTION FINISHED. P1,P2 REDEFINED
C  ALL OTHER SURFACE COEFFICIENTS ALSO REDEFINED

C  NOW: GENERAL SECOND ORDER EQUATION, RLB=1., FOR SURFACE I
C
2       CONTINUE
1     CONTINUE
C
      RETURN
990   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. SETFIT '
      WRITE (6,*) 'INCONSISTENCY IN IGNORABLE CO-ORDINATES DETECTED'
      WRITE (6,*) 'BETWEEN REQUESTING SURFACE I= ',I,' AND IE= ',IE
      WRITE (6,*) 'JUMLIM(I),LINFX,LINFY,LINFZ ',
     .             JUMLIM(I),LINFX,LINFY,LINFZ
      WRITE (6,*) 'EXIT CALLED '
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. SETFIT '
      WRITE (6,*) 'FIT OPTION FOR RLB(IE) = ',RLB(IE),' NOT FORESEEN'
      WRITE (6,*) 'REQUEST FROM SURFACE NO. ',I
      WRITE (6,*) 'IE = ',IE,' EXIT CALLED '
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. SETFIT '
      WRITE (6,*) 'THE VALID AREAS DO NOT INTERSECT'
      WRITE (6,*) 'I,IE ',I,IE
      CALL EXIT_OWN(1)
993   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. SETFIT '
      WRITE (6,*) 'STRAIGHT LINE NO I= ',I,' COLLAPSED TO A POINT'
      WRITE (6,*) 'SURFACE NO. I IS REDUNDANT. USE CH0 I/I OPTION '
      CALL EXIT_OWN(1)
      END
