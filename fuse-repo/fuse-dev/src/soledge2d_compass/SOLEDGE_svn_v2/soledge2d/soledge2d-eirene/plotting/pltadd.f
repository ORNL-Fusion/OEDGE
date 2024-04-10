C
C
      SUBROUTINE EIRENE_PLTADD (MANF,MEND)
C
C  PLOT ADDITIONAL SURFACES ISURF=MANF,MEND, PROJECTED INTO A PLANE
C  IF LZR,     PLOT IS DONE IN THIS ROUTINE
C  IF.NOT.LZR, PLOT DATA ARE COLLECTED IN MODULE CRECH
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CRECH
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLMSUR
      USE EIRMOD_CPLOT
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTEXT
      USE EIRMOD_CLGIN
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: MANF, MEND
      REAL(DP) :: G(4,2), R(4,2), AFF(3,3), AFFI(3,3)
      REAL(DP) :: EIRENE_ELLO, EIRENE_ELLU, EIRENE_PARA1, EIRENE_PARA2O,
     .          EIRENE_PARA2U,
     .          EIRENE_HYPP, EIRENE_HYPM, EIRENE_GERADY, EIRENE_GERADX
      REAL(DP) :: XTRAN, YTRAN, XI, ETA, XA, XE, XN, X1, X2,
     .          XANF, XEND, RAD, X13, Y13, X14, Y14, X23,
     .          Y23, X24, Y24, YANF, YEND, ARC05, XINI, XEN,
     .          ARC, YINI, XAN, XP, YP, X, Y, XS, YS, SQ, SQR,
     .          VR1, VR2, VR3, AMNPY, AMNPZ, AMXPY, AMXPZ, XMU,
     .          EIRENE_FMU, XP1, V1, V2, V3, XCHL, XCHR, YCHL, YCHR,
     .          DCHX, AMXPX, AMNPX, DCHY, YP1, TANA, SINAQ,
     .          COSAQ, ALF, DELS, DET, DETS, EIRENE_DETER,
     .          SCA, DH, EH, CK, DK, EK, BK, XP2, YP2, AK, AKEK2,
     .          S, DEL, AKEK, FK, DKCK, DKCK2
      INTEGER :: INN1, INN2, INN3, I, INN4, INN, IS, ICOLOR, J,
     .           ICUT, ISURF
      LOGICAL :: L1, L2, L3, L4, L5
      LOGICAL :: LBOX, LLISTE
      CHARACTER(1) :: CH1
      CHARACTER(2) :: CH2
      CHARACTER(3) :: CH3
      TYPE(PPOINT), POINTER :: CUR, SURFAN, SURFEN
 
      EXTERNAL EIRENE_ELLO,EIRENE_ELLU,EIRENE_PARA1,EIRENE_PARA2O,
     .         EIRENE_PARA2U,EIRENE_HYPP,EIRENE_HYPM,
     .         EIRENE_GERADY,EIRENE_GERADX
 
      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
C
      XCHL=CH2X0-CH2MX
      XCHR=CH2X0+CH2MX
      YCHL=CH2Y0-CH2MY
      YCHR=CH2Y0+CH2MY
      DCHX=2.*CH2MX
      DCHY=2.*CH2MY
      G(1,1)=XCHL
      G(1,2)=YCHL
      G(2,1)=XCHL
      G(2,2)=YCHR
      G(3,1)=XCHR
      G(3,2)=YCHL
      G(4,1)=XCHL
      G(4,2)=YCHL
      R(1,1)=0.
      R(1,2)=DCHY
      R(2,1)=DCHX
      R(2,2)=0.
      R(3,1)=0.
      R(3,2)=DCHY
      R(4,1)=DCHX
      R(4,2)=0.
C
C  PLOTTING PLANE IS X-Y, AT Z=CH2Z0
      DO 2000 ICUT=1,3
      IF (.NOT.PLCUT(ICUT)) GOTO 2000
      IF (ICUT.EQ.1) THEN
C  PLOT Y-Z AT X=CH2Z0. HENCE MAP:
C  Y --> Y
C  Z --> X
C  X --> Z
        IF (CH2Z0.NE.0.) CALL EIRENE_XSHADD(-CH2Z0,MANF,MEND)
        AFF=0.D0
        AFF(1,3)=1.
        AFF(2,2)=1.
        AFF(3,1)=-1.
        AFFI=0.D0
        AFFI(1,3)=-1.
        AFFI(2,2)=1.
        AFFI(3,1)=1.
        CALL EIRENE_ROTADD(AFF,AFFI,MANF,MEND)
      ELSEIF (ICUT.EQ.2) THEN
C  PLOT X-Z AT Y=CH2Z0. HENCE MAP:
C  X --> X
C  Z --> Y
C  Y --> -Z
        IF (CH2Z0.NE.0.) CALL EIRENE_YSHADD(-CH2Z0,MANF,MEND)
        AFF=0.D0
        AFF(1,1)=1.
        AFF(2,3)=1.
        AFF(3,2)=-1.
        AFFI=0.D0
        AFFI(1,1)=1.
        AFFI(2,3)=-1.
        AFFI(3,2)=1.
        CALL EIRENE_ROTADD(AFF,AFFI,MANF,MEND)
      ELSEIF (ICUT.EQ.3) THEN
C  PLOT X-Y AT Z=CH2Z0. HENCE MAP:
        IF (CH2Z0.NE.0.) CALL EIRENE_ZSHADD(-CH2Z0,MANF,MEND)
C  NO ROTATION NEEDED
      ENDIF
C
C
C   PLOT SURFACES
      ZPLT=0.
      INSTOR=0
C
      DO 200 ISURF=MANF,MEND
        IF (ILCOL(ISURF) == 666) GOTO 200
        J=ISURF
        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) INOSF=ISURF
C
C  PLOT SURFACE NO. ISURF
C
        IF (IGJUM0(ISURF).NE.0) GOTO 200
        IF (TRCPLT) THEN
          CALL EIRENE_LEER(1)
          WRITE (iunout,*) TXTSFL(ISURF)
        ENDIF
        LBOX=ILBOX(ISURF).GT.0.AND.LZR
        IF (LZR) THEN
          ICOLOR=ILCOL(ISURF)
          CALL GRNWPN(ICOLOR)
          IF (ILIIN(ISURF).LE.0) THEN
            CALL GRDSH(0.2,0.5,0.2)
          ELSE
            CALL GRDSH(1.,0.,1.)
          ENDIF
        ENDIF
C
C   PLOT SCHNITT MIT ALLGEMEINEM 3, 4 ODER 5 ECK
C
        IF (RLB(ISURF).GE.2.) THEN
          AMXPX=MAX(P1(1,ISURF),P2(1,ISURF))
          AMNPX=MIN(P1(1,ISURF),P2(1,ISURF))
          AMXPY=MAX(P1(2,ISURF),P2(2,ISURF))
          AMNPY=MIN(P1(2,ISURF),P2(2,ISURF))
          AMXPZ=MAX(P1(3,ISURF),P2(3,ISURF))
          AMNPZ=MIN(P1(3,ISURF),P2(3,ISURF))
          IF (RLB(ISURF).GE.3.) THEN
            AMXPX=MAX(AMXPX,P3(1,ISURF))
            AMNPX=MIN(AMNPX,P3(1,ISURF))
            AMXPY=MAX(AMXPY,P3(2,ISURF))
            AMNPY=MIN(AMNPY,P3(2,ISURF))
            AMXPZ=MAX(AMXPZ,P3(3,ISURF))
            AMNPZ=MIN(AMNPZ,P3(3,ISURF))
          ENDIF
          IF (RLB(ISURF).GE.4.) THEN
            AMXPX=MAX(AMXPX,P4(1,ISURF))
            AMNPX=MIN(AMNPX,P4(1,ISURF))
            AMXPY=MAX(AMXPY,P4(2,ISURF))
            AMNPY=MIN(AMNPY,P4(2,ISURF))
            AMXPZ=MAX(AMXPZ,P4(3,ISURF))
            AMNPZ=MIN(AMNPZ,P4(3,ISURF))
          ENDIF
          IF (RLB(ISURF).GE.5.) THEN
            AMXPX=MAX(AMXPX,P5(1,ISURF))
            AMNPX=MIN(AMNPX,P5(1,ISURF))
            AMXPY=MAX(AMXPY,P5(2,ISURF))
            AMNPY=MIN(AMNPY,P5(2,ISURF))
            AMXPZ=MAX(AMXPZ,P5(3,ISURF))
            AMNPZ=MIN(AMNPZ,P5(3,ISURF))
          ENDIF
          IF(AMXPX.LT.XCHL.OR.AMNPX.GT.XCHR.OR.
     .       AMXPY.LT.YCHL.OR.AMNPY.GT.YCHR.OR.
     .       AMXPZ.LT.ZPLT.OR.AMNPZ.GT.ZPLT) THEN
            IF (TRCPLT) WRITE (iunout,6660)
            GOTO 200
          ENDIF
          IF (TRCPLT) THEN
            WRITE (iunout,*) ' AMNPX,AMNPY,AMNPZ ',AMNPX,AMNPY,AMNPZ
            WRITE (iunout,*) ' AMXPX,AMXPY,AMXPZ ',AMXPX,AMXPY,AMXPZ
            WRITE (iunout,*) ' A0,A1,A2,A3 ',
     .                   A0LM(ISURF),A1LM(ISURF),A2LM(ISURF),A3LM(ISURF)
          ENDIF
c
          IF (ABS(A3LM(ISURF)).GE.1.-EPS10) THEN
C   N ECK PARALLEL ZU Z=ZPLT
            IF (TRCPLT) WRITE (iunout,*) 'N ECK PARALLEL ZU Z=ZPLT'
            IF (ABS(A0LM(ISURF)-ZPLT).GT.EPS10) GOTO 200
C
            L1=P1(1,ISURF).GE.XCHL.AND.P1(1,ISURF).LE.XCHR.AND.
     .         P1(2,ISURF).GE.YCHL.AND.P1(2,ISURF).LE.YCHR
            L2=P2(1,ISURF).GE.XCHL.AND.P2(1,ISURF).LE.XCHR.AND.
     .         P2(2,ISURF).GE.YCHL.AND.P2(2,ISURF).LE.YCHR
            L3=.FALSE.
            IF (RLB(ISURF).GE.3.) THEN
              L3=P3(1,ISURF).GE.XCHL.AND.P3(1,ISURF).LE.XCHR.AND.
     .           P3(2,ISURF).GE.YCHL.AND.P3(2,ISURF).LE.YCHR
            ENDIF
            L4=.FALSE.
            IF (RLB(ISURF).GE.4.) THEN
              L4=P4(1,ISURF).GE.XCHL.AND.P4(1,ISURF).LE.XCHR.AND.
     .           P4(2,ISURF).GE.YCHL.AND.P4(2,ISURF).LE.YCHR
            ENDIF
            L5=.FALSE.
            IF (RLB(ISURF).GE.5.) THEN
              L5=P5(1,ISURF).GE.XCHL.AND.P5(1,ISURF).LE.XCHR.AND.
     .           P5(2,ISURF).GE.YCHL.AND.P5(2,ISURF).LE.YCHR
            ENDIF
            IF (TRCPLT) WRITE (iunout,*) 'L1,L2,L3,L4,L5 ',
     .                                    L1,L2,L3,L4,L5
C   PLOT P1, P2
            IF (L1.AND.L2) THEN
              IF (LZR) THEN
                CALL GRJMP (REAL(P1(1,ISURF),KIND(1.E0)),
     .                      REAL(P1(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P2(1,ISURF),KIND(1.E0)),
     .                      REAL(P2(2,ISURF),KIND(1.E0)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL EIRENE_STCOOR (P1(1,ISURF),P1(2,ISURF),0)
                CALL EIRENE_STCOOR (P2(1,ISURF),P2(2,ISURF),1)
              ENDIF
            ELSE
              CALL EIRENE_PSIDE
     .  (G,R,P1(1,ISURF),P1(2,ISURF),P2(1,ISURF),
     .                    P2(2,ISURF),L1,L2,EPS10)
            ENDIF
C   PLOT P1, P3
            IF (L1.AND.L3) THEN
              IF (LZR) THEN
                CALL GRJMP (REAL(P1(1,ISURF),KIND(1.E0)),
     .                      REAL(P1(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P3(1,ISURF),KIND(1.E0)),
     .                      REAL(P3(2,ISURF),KIND(1.E0)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL EIRENE_STCOOR (P1(1,ISURF),P1(2,ISURF),0)
                CALL EIRENE_STCOOR (P3(1,ISURF),P3(2,ISURF),1)
              ENDIF
            ELSE
              CALL EIRENE_PSIDE
     .  (G,R,P1(1,ISURF),P1(2,ISURF),P3(1,ISURF),
     .                    P3(2,ISURF),L1,L3,EPS10)
            ENDIF
C   PLOT P2, P3
            IF (L2.AND.L3.AND..NOT.L4) THEN
              IF (LZR) THEN
                CALL GRJMP (REAL(P2(1,ISURF),KIND(1.E0)),
     .                      REAL(P2(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P3(1,ISURF),KIND(1.E0)),
     .                      REAL(P3(2,ISURF),KIND(1.E0)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL EIRENE_STCOOR (P2(1,ISURF),P2(2,ISURF),0)
                CALL EIRENE_STCOOR (P3(1,ISURF),P3(2,ISURF),1)
              ENDIF
            ELSE
              CALL EIRENE_PSIDE
     .  (G,R,P2(1,ISURF),P2(2,ISURF),P3(1,ISURF),
     .                    P3(2,ISURF),L2,L3,EPS10)
            ENDIF
C   PLOT P4, P2
            IF (L4.AND.L2) THEN
              IF (LZR) THEN
                CALL GRJMP (REAL(P4(1,ISURF),KIND(1.E0)),
     .                      REAL(P4(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P2(1,ISURF),KIND(1.E0)),
     .                      REAL(P2(2,ISURF),KIND(1.E0)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL EIRENE_STCOOR (P4(1,ISURF),P4(2,ISURF),0)
                CALL EIRENE_STCOOR (P2(1,ISURF),P2(2,ISURF),1)
              ENDIF
            ELSE
              CALL EIRENE_PSIDE
     .  (G,R,P4(1,ISURF),P4(2,ISURF),P2(1,ISURF),
     .                    P2(2,ISURF),L4,L2,EPS10)
            ENDIF
C   PLOT P4, P3
            IF (L4.AND.L3.AND..NOT.L5) THEN
              IF (LZR) THEN
                CALL GRJMP (REAL(P4(1,ISURF),KIND(1.E0)),
     .                      REAL(P4(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P3(1,ISURF),KIND(1.E0)),
     .                      REAL(P3(2,ISURF),KIND(1.E0)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL EIRENE_STCOOR (P4(1,ISURF),P4(2,ISURF),0)
                CALL EIRENE_STCOOR (P3(1,ISURF),P3(2,ISURF),1)
              ENDIF
            ELSE
              CALL EIRENE_PSIDE
     .  (G,R,P4(1,ISURF),P4(2,ISURF),P3(1,ISURF),
     .                    P3(2,ISURF),L4,L3,EPS10)
            ENDIF
C   PLOT P4, P5
            IF (L4.AND.L5) THEN
              IF (LZR) THEN
                CALL GRJMP (REAL(P4(1,ISURF),KIND(1.E0)),
     .                      REAL(P4(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P5(1,ISURF),KIND(1.E0)),
     .                      REAL(P5(2,ISURF),KIND(1.E0)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL EIRENE_STCOOR (P4(1,ISURF),P4(2,ISURF),0)
                CALL EIRENE_STCOOR (P5(1,ISURF),P5(2,ISURF),1)
              ENDIF
            ELSE
              CALL EIRENE_PSIDE
     .  (G,R,P4(1,ISURF),P4(2,ISURF),P5(1,ISURF),
     .                    P5(2,ISURF),L4,L5,EPS10)
            ENDIF
C   PLOT P5, P3
            IF (L5.AND.L3) THEN
              IF (LZR) THEN
                CALL GRJMP (REAL(P5(1,ISURF),KIND(1.E0)),
     .                      REAL(P5(2,ISURF),KIND(1.E0)))
                CALL GRDRW (REAL(P3(1,ISURF),KIND(1.E0)),
     .                      REAL(P3(2,ISURF),KIND(1.E0)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL EIRENE_STCOOR (P5(1,ISURF),P5(2,ISURF),0)
                CALL EIRENE_STCOOR (P3(1,ISURF),P3(2,ISURF),1)
              ENDIF
            ELSE
              CALL EIRENE_PSIDE
     .  (G,R,P5(1,ISURF),P5(2,ISURF),P3(1,ISURF),
     .                    P3(2,ISURF),L5,L3,EPS10)
            ENDIF
C   BERECHNE SCHNITTGERADE
          ELSE
            VR1=-A2LM(ISURF)
            VR2=A1LM(ISURF)
            VR3=0.
            V3=ZPLT
            IF (TRCPLT) WRITE (iunout,*) ' VR1,VR2,VR3 ',VR1,VR2,VR3
            IF (ABS(A1LM(ISURF)).GT.ABS(A2LM(ISURF))) THEN
              V1=(-A0LM(ISURF)-A3LM(ISURF)*V3)/A1LM(ISURF)
              V2=0.
            ELSE
              V1=0.
              V2=(-A0LM(ISURF)-A3LM(ISURF)*V3)/A2LM(ISURF)
            ENDIF
            IS=0
C
            XMU=EIRENE_FMU(V1,V2,V3,VR1,VR2,VR3,
     .              P1(1,ISURF),P1(2,ISURF),P1(3,ISURF),
     .              P2(1,ISURF)-P1(1,ISURF),P2(2,ISURF)-P1(2,ISURF),
     .              P2(3,ISURF)-P1(3,ISURF),EPS10)
            IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
              IS=1
              XP1=P1(1,ISURF)+XMU*(P2(1,ISURF)-P1(1,ISURF))
              YP1=P1(2,ISURF)+XMU*(P2(2,ISURF)-P1(2,ISURF))
            ENDIF
            IF (RLB(ISURF).GE.3.) THEN
              XMU=EIRENE_FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P1(1,ISURF),P1(2,ISURF),P1(3,ISURF),
     .                P3(1,ISURF)-P1(1,ISURF),P3(2,ISURF)-P1(2,ISURF),
     .                P3(3,ISURF)-P1(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P1(1,ISURF)+XMU*(P3(1,ISURF)-P1(1,ISURF))
                  YP1=P1(2,ISURF)+XMU*(P3(2,ISURF)-P1(2,ISURF))
                ELSE
                  XP2=P1(1,ISURF)+XMU*(P3(1,ISURF)-P1(1,ISURF))
                  YP2=P1(2,ISURF)+XMU*(P3(2,ISURF)-P1(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
            ENDIF
            IF (RLB(ISURF).GE.4.) THEN
              XMU=EIRENE_FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P4(1,ISURF),P4(2,ISURF),P4(3,ISURF),
     .                P2(1,ISURF)-P4(1,ISURF),P2(2,ISURF)-P4(2,ISURF),
     .                P2(3,ISURF)-P4(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P4(1,ISURF)+XMU*(P2(1,ISURF)-P4(1,ISURF))
                  YP1=P4(2,ISURF)+XMU*(P2(2,ISURF)-P4(2,ISURF))
                ELSE
                  XP2=P4(1,ISURF)+XMU*(P2(1,ISURF)-P4(1,ISURF))
                  YP2=P4(2,ISURF)+XMU*(P2(2,ISURF)-P4(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
              XMU=EIRENE_FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P4(1,ISURF),P4(2,ISURF),P4(3,ISURF),
     .                P3(1,ISURF)-P4(1,ISURF),P3(2,ISURF)-P4(2,ISURF),
     .                P3(3,ISURF)-P4(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P4(1,ISURF)+XMU*(P3(1,ISURF)-P4(1,ISURF))
                  YP1=P4(2,ISURF)+XMU*(P3(2,ISURF)-P4(2,ISURF))
                ELSE
                  XP2=P4(1,ISURF)+XMU*(P3(1,ISURF)-P4(1,ISURF))
                  YP2=P4(2,ISURF)+XMU*(P3(2,ISURF)-P4(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
            ENDIF
            IF (RLB(ISURF).GE.5) THEN
              XMU=EIRENE_FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P4(1,ISURF),P4(2,ISURF),P4(3,ISURF),
     .                P5(1,ISURF)-P4(1,ISURF),P5(2,ISURF)-P4(2,ISURF),
     .                P5(3,ISURF)-P4(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P4(1,ISURF)+XMU*(P5(1,ISURF)-P4(1,ISURF))
                  YP1=P4(2,ISURF)+XMU*(P5(2,ISURF)-P4(2,ISURF))
                ELSE
                  XP2=P4(1,ISURF)+XMU*(P5(1,ISURF)-P4(1,ISURF))
                  YP2=P4(2,ISURF)+XMU*(P5(2,ISURF)-P4(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
              XMU=EIRENE_FMU(V1,V2,V3,VR1,VR2,VR3,
     .                P5(1,ISURF),P5(2,ISURF),P5(3,ISURF),
     .                P3(1,ISURF)-P5(1,ISURF),P3(2,ISURF)-P5(2,ISURF),
     .                P3(3,ISURF)-P5(3,ISURF),EPS10)
              IF (XMU.GT.-EPS10.AND.XMU.LT.1.+EPS10) THEN
                IS=IS+1
                IF (IS.EQ.1) THEN
                  XP1=P5(1,ISURF)+XMU*(P3(1,ISURF)-P5(1,ISURF))
                  YP1=P5(2,ISURF)+XMU*(P3(2,ISURF)-P5(2,ISURF))
                ELSE
                  XP2=P5(1,ISURF)+XMU*(P3(1,ISURF)-P5(1,ISURF))
                  YP2=P5(2,ISURF)+XMU*(P3(2,ISURF)-P5(2,ISURF))
                  GOTO 90
                ENDIF
              ENDIF
            ENDIF
C
            IF (IS.EQ.0) THEN
              WRITE (iunout,*) ' SCHNITTGERADE AUSSERHALB DES VIERECKS'
            ELSEIF (IS.EQ.1.AND.
     .              XP1.GE.XCHL.AND.XP1.LE.XCHR.AND.YP1.GE.YCHL.AND.
     .              YP1.LE.YCHR) THEN
              IF (LZR) THEN
                CALL GRSPTS (30)
                CALL GRJMP
     .  (REAL(XP1,KIND(1.E0)),REAL(YP1,KIND(1.E0)))
                CALL GRDRW
     .  (REAL(XP1,KIND(1.E0)),REAL(YP1,KIND(1.E0)))
                CALL GRSPTS(16)
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL EIRENE_STCOOR
     .  (XP1,YP1,0)
            ENDIF
            GOTO 200
C
90          CONTINUE
            L1=XP1.GE.XCHL.AND.XP1.LE.XCHR.AND.
     .         YP1.GE.YCHL.AND.YP1.LE.YCHR
            L2=XP2.GE.XCHL.AND.XP2.LE.XCHR.AND.
     .         YP2.GE.YCHL.AND.YP2.LE.YCHR
            IF (L1.AND.L2) THEN
              IF (LZR) THEN
                CALL GRJMP
     .  (REAL(XP1,KIND(1.E0)),REAL(YP1,KIND(1.E0)))
                CALL GRDRW
     .  (REAL(XP2,KIND(1.E0)),REAL(YP2,KIND(1.E0)))
              ENDIF
              IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                CALL EIRENE_STCOOR (XP1,YP1,0)
                CALL EIRENE_STCOOR (XP2,YP2,1)
              ENDIF
            ELSE
              CALL EIRENE_PSIDE (G,R,XP1,YP1,XP2,YP2,L1,L2,EPS10)
            ENDIF
          ENDIF
          GOTO 200
        ENDIF
C
C   ALLE ANDEREN FALLE: RLB.LT.2.
C
        IF (RLB(ISURF).EQ.1.) THEN
          IF (ZLIMS1(1,ISURF).GT.ZPLT.OR.ZLIMS2(1,ISURF).LT.ZPLT) THEN
            IF (TRCPLT) THEN
              WRITE (iunout,*) 'ZLIMS1,ZPLT,ZLIMS2 ',
     .                     ZLIMS1(1,ISURF),ZPLT,ZLIMS2(1,ISURF)
            ENDIF
            GOTO 200
          ENDIF
        ENDIF
        XM=0.
        YM=0.
        AK=A4LM(ISURF)
        BK=A7LM(ISURF)
        CK=A5LM(ISURF)
        DK=A1LM(ISURF)+A8LM(ISURF)*ZPLT
        EK=A2LM(ISURF)+A9LM(ISURF)*ZPLT
        FK=A0LM(ISURF)+(A3LM(ISURF)+A6LM(ISURF)*ZPLT)*ZPLT
        IF (TRCPLT) THEN
          WRITE (iunout,*) 'COEFFICIENTS OF CURVE-EQUATION'
          WRITE (iunout,*) 'AK*X**2+BK*XY+CK*Y**2+DK*X+EK*Y+FK=0.    '
          WRITE (iunout,6662)  AK,BK,CK,DK,EK,FK
        ENDIF
C
C   DIE KURVENGLEICHUNG LAUTET
C   AK*X**2+BK*XY+CK*Y**2+DK*X+EK*Y+FK=0.    '
C   TRANSFORMIERE BK AUF 0., ABER ERHALTE INVARIANTEN
C
        DKCK=DK*CK
        DKCK2=DKCK+DKCK
        AKEK=AK*EK
        AKEK2=AKEK+AKEK
C  ALTE INVARIANTEN:
        S=AK+CK
        DEL=AK*CK-BK*BK*0.25
        DELS=DEL
        IF (ABS(DEL).LE.EPS10) DEL=0.
        DET=FK*DEL+(DK*(BK*EK-DKCK2)-EK*(AKEK2-DK*BK))*0.125
        DETS=DET
        IF (ABS(DET).LE.EPS10) DET=0.
        IF (TRCPLT) THEN
          WRITE (iunout,*) 'S,DELS,DETS,DEL,DET '
          WRITE (iunout,*) S,DELS,DETS,DEL,DET
        ENDIF
C
        ALF=0.
        IF (BK.EQ.0.) THEN
          COSA=1.
          SINA=0.
          TANA=0.
          A=AK
          C=CK
          D=DK
          E=EK
          F=FK
        ELSE
          IF (AK.NE.CK) ALF=ATAN(BK/(AK-CK))*0.5
          IF (AK.EQ.CK) ALF=45.*PIA/180.
          SINA=SIN(ALF)
          IF (SIGN(1._DP,SINA).NE.SIGN(1._DP,BK)) THEN
            ALF=ALF+PIA
            SINA=SIN(ALF)
          ENDIF
          COSA=COS(ALF)
          TANA=TAN(ALF)
          SINAQ=SINA*SINA
          COSAQ=1.-SINAQ
          SCA=SINA*COSA
          A=AK*COSAQ+BK*SCA+CK*SINAQ
          IF (ABS(A).LT.EPS10) A=0.
          C=S-A
          IF (ABS(C).LT.EPS10) C=0.
          F=FK
          IF (ABS(F).LT.EPS10) F=0.
          D=DK*COSA+EK*SINA
          IF (ABS(D).LT.EPS10) D=0.
          E=EK*COSA-DK*SINA
          IF (ABS(E).LT.EPS10) E=0.
          DH=D*0.5
          EH=E*0.5
          DET=EIRENE_DETER(A,0._DP,DH,0._DP,C,EH,DH,EH,F)
          IF (ABS(DET).LE.EPS10) DET=0.
        ENDIF
C
        IF (TRCPLT) THEN
          WRITE (iunout,*)
     .      'AFTER TRANSFORMATION: B=0., NEW COEFFICIENTS'
          WRITE (iunout,6661) A,C,D,E,F
          WRITE (iunout,*) 'DETERMINANT DET= ',DET
        ENDIF
C
        IF (RLB(ISURF).EQ.1.) THEN
          XL1=MAX(XCHL,XLIMS1(1,ISURF))
          XL1=MIN(XCHR,XL1)
          XL2=MAX(XCHL,XLIMS2(1,ISURF))
          XL2=MIN(XCHR,XL2)
          YL1=MAX(YCHL,YLIMS1(1,ISURF))
          YL1=MIN(YCHR,YL1)
          YL2=MAX(YCHL,YLIMS2(1,ISURF))
          YL2=MIN(YCHR,YL2)
        ELSE
          XL1=XCHL
          XL2=XCHR
          YL1=YCHL
          YL2=YCHR
        ENDIF
C
        X13=YL1*SINA+XL1*COSA
        X14=YL2*SINA+XL1*COSA
        X23=YL1*SINA+XL2*COSA
        X24=YL2*SINA+XL2*COSA
C
        XANF=MIN(X13,X14,X23,X24)
        XEND=MAX(X13,X14,X23,X24)
C
        Y13=-XL1*SINA+YL1*COSA
        Y14=-XL2*SINA+YL1*COSA
        Y23=-XL1*SINA+YL2*COSA
        Y24=-XL2*SINA+YL2*COSA
C
        YANF=MIN(Y13,Y14,Y23,Y24)
        YEND=MAX(Y13,Y14,Y23,Y24)
C
        IF (RLB(ISURF).EQ.1.5) THEN
          XC1=XLIMS1(1,ISURF)
          XC2=XLIMS2(1,ISURF)
          YC1=YLIMS1(1,ISURF)
          YC2=YLIMS2(1,ISURF)
        ENDIF
C
        IF (RLB(ISURF).LT.0.) THEN
           DO 401 I=1,ILIN(ISURF)
              ALIN(I)=ALIMS(I,ISURF)
              XLIN(I)=XLIMS(I,ISURF)
              YLIN(I)=YLIMS(I,ISURF)
              ZLIN(I)=ZLIMS(I,ISURF)
401        CONTINUE
           DO 402 I=1,ISCN(ISURF)
              A0S(I)=ALIMS0(I,ISURF)
              A1S(I)=XLIMS1(I,ISURF)
              A2S(I)=YLIMS1(I,ISURF)
              A3S(I)=ZLIMS1(I,ISURF)
              A4S(I)=XLIMS2(I,ISURF)
              A5S(I)=YLIMS2(I,ISURF)
              A6S(I)=ZLIMS2(I,ISURF)
              A7S(I)=XLIMS3(I,ISURF)
              A8S(I)=YLIMS3(I,ISURF)
              A9S(I)=ZLIMS3(I,ISURF)
402        CONTINUE
           MLIN=ILIN(ISURF)
           MSCN=ISCN(ISURF)
        ENDIF
        IF (TRCPLT) THEN
          WRITE (iunout,*) 'PLOT REGION:'
          WRITE (iunout,*) 'XL1,XL2,YL1,YL2 ',XL1,XL2,YL1,YL2
          WRITE (iunout,*) 'PLOT REGION AFTER TRANSFORMATION'
          WRITE (iunout,*) 'XANF,XEND,YANF,YEND ',XANF,XEND,YANF,YEND
        ENDIF
C
C
C     FALLUNTERSCHEIDUNGEN
C
        IF (ABS(DEL).LE.EPS10) GOTO 1000
        XN=-DET/DEL
        IF (ABS(XN).LE.EPS10) XN=0.
        IF (XN) 405,500,400
405     IF (DEL.LT.-EPS10) GOTO 450
        IF (A.LT.-EPS10.AND.C.LT.-EPS10) GOTO 410
        GOTO 407
400     IF (DEL.LT.-EPS10) GOTO 450
        IF (A.GT.EPS10.AND.C.GT.EPS10) GOTO 410
C
C     KEINE REELLE LOESUNG
C
407     IF (TRCPLT) WRITE (iunout,6664)
        GOTO 200
C
C     ELLIPSE
C
410     AHALB=SQRT(ABS(XN/A))
        BHALB=SQRT(ABS(XN/C))
        XM=-D/(2.*A)
        YM=-E/(2.*C)
        IF (RLB(ISURF).EQ.1.5) THEN
           XA=-AHALB
           XE=AHALB
        ELSE
C   PLOTGEBIET ENTHAELT GEDREHTE BEGRENZUNGSBOX
C   BEGRENZUNGSBOX WIRD ALS PLOTGRENZE NACH RUECKTRANSFORMATION GENUTZT
           XA=MAX(XM-AHALB,XANF+1.E-8_DP)-XM
           XE=MIN(XM+AHALB,XEND-1.E-8_DP)-XM
        ENDIF
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_ELLO,XA,XE,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'ELLO'
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_ELLU,XA,XE,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'ELLU'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6665)
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6665)
        GOTO 200
C
C     HYPERBEL
C
450     CONTINUE
        RAD=MAX(0._DP,XN/A)
        X1=-D/A/2.+SQRT(RAD)
        X2=-D/A/2.-SQRT(RAD)
        INN1=0
        INN2=0
        INN3=0
        INN4=0
        IF (X1.GT.XL2) GOTO  451
        XA=MAX(X1,XANF)+1.D-6
        XE=XEND
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_HYPP,XA,XE,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'HYPP'
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_HYPM,XA,XE,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'HYPM'
        IF (INN1+INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6666)
451     IF (X2.LT.XL1) GOTO  452
        XA=XANF
        XE=MIN(X2,XEND)-1.D-6
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_HYPP,XA,XE,INN3,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'HYPP'
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_HYPM,XA,XE,INN4,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'HYPM'
452     IF (INN3+INN4.EQ.0.AND.TRCPLT) WRITE (iunout,6666)
        GOTO 200
C
C
C     EIN PUNKT
C
500     CONTINUE
        IF (DEL.LT.-EPS10) GOTO 510
        X=-D/(2.*A)
        Y=-E/(2.*C)
        XP=XTRAN(X,Y)
        YP=YTRAN(X,Y)
        IF (TRCPLT) WRITE (iunout,*) 'POINT'
        IF (XP.LT.XL1.OR.XP.GT.XL2.OR.YP.LT.YL1.OR.YP.GT.YL2) THEN
          IF (TRCPLT) WRITE (iunout,6667)
          GOTO 200
        ENDIF
        IF (LZR) THEN
          CALL GRCHRC (0.1,0.0,16)
          CALL GRJMPS (REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)),3)
          CALL GRCHRC (0.3,0.0,16)
        ENDIF
        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL EIRENE_STCOOR (XP,YP,0)
        GOTO 200
C
C     ZWEI SICH SCHNEIDENDE GERADEN
C
510     CONTINUE
        A0=-SQRT(ABS(A/C))
        A1=-SQRT(ABS(A/C))*D/(2.*A)-E/(2.*C)
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_GERADY,XANF,XEND,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
        A0=SQRT(ABS(A/C))
        A1=SQRT(ABS(A/C))*D/(2.*A)-E/(2.*C)
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_GERADY,XANF,XEND,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
        GOTO 200
C
1000    CONTINUE
        IF (ABS(A).GT.EPS10) GOTO 1300
        IF (ABS(C).GT.EPS10) GOTO 1400
        IF (ABS(E).GT.EPS10.AND.ABS(E).GE.ABS(D)) GOTO 1200
        IF (ABS(D).GT.EPS10.AND.ABS(D).GE.ABS(E)) GOTO 1100
        IF (ABS(F).GT.EPS10) GOTO 1002
        IF (TRCPLT) THEN
        WRITE (iunout,*) 'GLEICHUNG DER FORM 0.=0.'
        ENDIF
        GOTO 200
1002    IF (TRCPLT) WRITE (iunout,1001) F
1001    FORMAT (//1X,'KEIN PLOTT, DENN GLEICHUNG DER FORM F=',
     .            1PE12.4,' = 0')
        GOTO 200
C
C     GERADE;  D*X + E*Y + F = 0, D UNGLEICH 0., E=0. MOEGLICH
C
1100    CONTINUE
        A0=-E/D
        A1=-F/D
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_GERADX,YANF,YEND,INN,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADX'
        IF (INN.EQ.0.AND.TRCPLT) WRITE (iunout,6669)
        GOTO 200
C
C     GERADE; D*X + E*Y + F = 0, E UNGLEICH 0., D=0. MOEGLICH
C
1200    CONTINUE
        A0=-D/E
        A1=-F/E
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_GERADY,XANF,XEND,INN,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN.EQ.0.AND.TRCPLT) WRITE (iunout,6669)
        GOTO 200
C
C     PARABEL; A*X**2 + D*X + E*Y + F = 0
C
1300    CONTINUE
        IF (ABS(DET).LE.EPS10) GOTO 1310
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_PARA1,XANF,XEND,INN,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'PARA1'
        IF (INN.EQ.0.AND.TRCPLT) WRITE (iunout,6668)
        GOTO 200
C
C     PAAR PARALLELER GERADEN; A*X**2 + D*X + F = 0
C
1310    SQ=(D*D/(4.*A)-F)/A
        SQR=0.
        IF (SQ.GT.0.) SQR=SQRT(SQ)
        IF (SQ.LT.0.) THEN
          IF (TRCPLT) WRITE (iunout,6664)
          GOTO 200
        ENDIF
        A0=0.
        A1=-D/(2.*A)+SQR
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_GERADX,YANF,YEND,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADX'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
        A1=-D/(2.*A)-SQR
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_GERADX,YANF,YEND,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADX'
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
        GOTO 200
C
C     PARABEL; C*Y**2 + D*X + E*Y + F = 0
C
1400    CONTINUE
        IF (ABS(DET).LE.EPS10) GOTO 1500
        XS=(E*E/(4.*C)-F)/D
        YS=EIRENE_PARA2O(XS+1.)
        IF (YS.LT.1.D50) THEN
           XAN=MIN(MAX(XANF,XS),XEND)
           XEN=XEND
        ELSE
           XEN=MAX(MIN(XEND,XS),XANF)
           XAN=XANF
        ENDIF
        IF (ABS(XEN-XAN).LT.EPS10) THEN
           IF (TRCPLT) WRITE (iunout,6668)
           GOTO 200
        ENDIF
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_PARA2O,XAN,XEN,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'PARA2O'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6671)
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_PARA2U,XAN,XEN,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'PARA2U'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6671)
        GOTO 200
C
C     PAAR PARALLELER GERADEN; C*Y**2 + E*Y + F = 0
C
1500    CONTINUE
        A0=0.
        SQ=(E*E/(4.*C)-F)/C
        SQR=0
        IF (SQ.GT.0.) SQR=SQRT(SQ)
        IF (SQ.LT.0.) THEN
           IF (TRCPLT) WRITE (iunout,6664)
           GOTO 200
        ENDIF
        A1=-E/(2.*C)+SQR
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_GERADY,XANF,XEND,INN1,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN1.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
        A1=-E/(2.*C)-SQR
        CALL EIRENE_PLTIN (RLB(ISURF),EIRENE_GERADY,XANF,XEND,INN2,LBOX)
        IF (TRCPLT) WRITE (iunout,*) 'GERADY'
        IF (INN2.EQ.0.AND.TRCPLT) WRITE (iunout,6670)
C
C  SURFACE NO ISURF DONE
C
200   CONTINUE
C
C  END OF DO LOOP OVER SURFACE NUMBERS ISURF
C
      IF (LZR) THEN
        CALL GRDSH(1.,0.,1.)
        CALL GRNWPN(1)
      ENDIF
C
      IF (ICUT.EQ.1) THEN
        CALL EIRENE_ROTADD(AFFI,AFF,MANF,MEND)
        IF (CH2Z0.NE.0.) CALL EIRENE_XSHADD(CH2Z0,MANF,MEND)
      ELSEIF (ICUT.EQ.2) THEN
        CALL EIRENE_ROTADD(AFFI,AFF,MANF,MEND)
        IF (CH2Z0.NE.0.) CALL EIRENE_YSHADD(CH2Z0,MANF,MEND)
      ELSEIF (ICUT.EQ.3) THEN
        IF (CH2Z0.NE.0.) CALL EIRENE_ZSHADD(CH2Z0,MANF,MEND)
C  NO INVERS ROTATION NEEDED
      ENDIF
C
2000  CONTINUE
C
      IF (TRCPLT) WRITE (iunout,*) 'INSTOR= ',INSTOR
C
      IF (PLNUMS.OR.PLARR) THEN
        CUR => FIRST_POINT
        DO WHILE (ASSOCIATED(CUR))
          IF (CUR%NPL2D == 0) SURFAN => CUR
          ARC=0.
          lliste = ASSOCIATED(CUR%NXTPNT)
          if (lliste) lliste = lliste .and. (CUR%NXTPNT%NPL2D == 1)
          DO WHILE (lliste)
            ARC=ARC+SQRT((CUR%XPL2D - CUR%NXTPNT%XPL2D)**2 +
     .                   (CUR%YPL2D - CUR%NXTPNT%YPL2D)**2)
            CUR => CUR%NXTPNT
            lliste = lliste .and. ASSOCIATED(CUR%NXTPNT)
            if (lliste) lliste = lliste .and. (CUR%NXTPNT%NPL2D == 1)
          END DO
          SURFEN => CUR
          ARC05=ARC*0.5
          ARC=0.
          CUR => SURFAN
          DO WHILE (ARC < ARC05)
            ARC=ARC+SQRT((CUR%XPL2D - CUR%NXTPNT%XPL2D)**2 +
     .                   (CUR%YPL2D - CUR%NXTPNT%YPL2D)**2)
            IF (ARC < ARC05) CUR => CUR%NXTPNT
          END DO
C  POINT BETWEEN CUR AND CUR%NXTPNT
          XINI=(CUR%XPL2D+CUR%NXTPNT%XPL2D)*0.5
          YINI=(CUR%YPL2D+CUR%NXTPNT%YPL2D)*0.5
C
C  PLOT ARROWS: SURFACE NORMAL
          IF (PLARR) THEN
c           xlst=
c           ylst=
c           alen=
c           awid=
c           icode=
c           call grarrw(REAL(XINI,KIND(1.E0)),REAL(YINI,KIND(1.E0)),xlst,ylst,
c    .                  REAL(alen,KIND(1.E0)),REAL(awid,KIND(1.E0)),icode)
c
          ENDIF
C  PLOT SURFACE NUMBERS
          IF (PLNUMS) THEN
            IF (CUR%NUMSUR.LT.10) THEN
              WRITE (CH1,'(I1)') CUR%NUMSUR
              CALL GRTXT
     .  (REAL(XINI,KIND(1.E0)),REAL(YINI,KIND(1.E0)),
     .                    1,CH1)
            ELSEIF (CUR%NUMSUR.LT.100) THEN
              WRITE (CH2,'(I2)') CUR%NUMSUR
              CALL GRTXT
     .  (REAL(XINI,KIND(1.E0)),REAL(YINI,KIND(1.E0)),
     .                    2,CH2)
            ELSEIF (CUR%NUMSUR.LT.1000) THEN
              WRITE (CH3,'(I3)') CUR%NUMSUR
              CALL GRTXT
     .  (REAL(XINI,KIND(1.E0)),REAL(YINI,KIND(1.E0)),
     .                    3,CH3)
            ENDIF
          ENDIF
          CUR => SURFEN%NXTPNT
        END DO
      ENDIF
C

      call eirene_vtkout

      RETURN
6660  FORMAT (//1X,' CLOSED POLYGON OUTSIDE PLOTREGION')
6661  FORMAT (/1X,' A,C,D,E,F',/1X,1P,5E12.4)
6662  FORMAT (/1X,' AK,BK,CK,DK,EK,FK',/1X,1P,6E12.4)
6664  FORMAT (//1X,'NO REELL SOLUTION')
6665  FORMAT (//1X,'HALF ELLIPSE OUTSIDE PLOTREGION')
6666  FORMAT (//1X,'HYPERBEL OUTSIDE PLOTREGION')
6667  FORMAT (//1X,'SINGLE POINT OUTSIDE PLOTREGION')
6668  FORMAT (//1X,'PARABEL OUTSIDE PLOTREGION')
6669  FORMAT (//1X,'STRAIGHT LINE OUTSIDE PLOTREGION')
6670  FORMAT (//1X,'ONE OF THE TWO STRAIGHT LINES OUTSIDE PLOTREGION')
6671  FORMAT (//1X,'HALF PARABEL OUTSIDE PLOTREGION')
      END
