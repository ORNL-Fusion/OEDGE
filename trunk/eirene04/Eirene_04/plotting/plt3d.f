C
C
      SUBROUTINE PLT3D (XR,YR,FAKX,FAKY,ITH,ABSMIN,ABSMAX,ORDMIN,ORDMAX)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CRECH
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CPL3D
      USE CPLOT
      USE CINIT
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CTRCEI
      USE CGEOM
      USE COMSOU
      USE CTEXT
      USE CLGIN

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: XR, YR, FAKX, FAKY,
     .                       ABSMIN, ABSMAX, ORDMIN, ORDMAX
      INTEGER, INTENT(OUT) :: ITH

      INTEGER, PARAMETER :: NPLY=501
      REAL(DP) :: AL(10), AR(10), XP(NPLY), YP(NPLY), ZPLOT(N3RD+NTOR),
     .          XSAVE(NPLY,N3RD+NTOR), YSAVE(NPLY,N3RD+NTOR),
     .          PHIAN(9),PHIEN(9)
      REAL(DP) :: XX(101),YY(101)
      REAL(DP) :: DM, RS, Y, TR, EP, EL, PHI, TA, TD, F1B, F2B, F3B,
     .          TB, RR, X, Z, CH2MXS, CH2MYS,
     .          CH2X0S, CH2Y0S, XCENT, YCENT, XNULL, YNULL, DELX, DELY,
     .          XH, YH, XMA3D, YMA3D, ZMA3D, XMI3D, YMI3D, ZMI3D, XN,
     .          YN, ZX0B, ZY0B, ZZ0B, PHA, PHE, CXB, B1B, B2B, B3B,
     .          CYB, CZB, F0B, B0B, RZYLB, CX, CY, CZ, ZX0, ZY0, ZZ0,
     .          RZYL, T1, T2, F0, F1, F2, F3
      REAL(SP) :: XPS(NPLY),YPS(NPLY)
      INTEGER:: ILT(N3RD+NTOR)
      INTEGER :: IR, IBR, IST, IS, NR, I1, IZ, NP, II, K, ID, NA, ISTP,
     .           IA, IAN, IEN, KIN, ISSTD, IBA, NJZ, J, JJ, IPZ,
     .           IP, I, NINNE, NZAD, NIN, MERK2, IPR, IB, MERK,
     .           IJZ, JP, IPRT
      LOGICAL :: PLABLE(NLIM), LPERID(NLIM), LSYMET(NLIM),
     .           LERR1, LERR2, LSAVE, PLT1, PLT2, PLT3
      TYPE(PPOINT), POINTER :: CUR
C
C  DEFAULT TOROIDAL GRID PLOT OPTION
C  IN CASE THAT NO TOROIDAL GRID IS DEFINED
C
      IF (.NOT.NLTOR.AND..NOT.NLTRA) THEN
        IPLTS(3)=1
        IPLAS(3,1)=1
        IPLES(3,1)=2
      ENDIF
C
      IF (NLTRA) THEN
        ITH=ANGLE3
        IF (ITH.LE.0.OR.ITH.GE.NTTRA) ITH=1
        WIN=ZZONE(ITH)
        RMT=RMTOR
        WINJ=WIN
      ELSEIF (NLTRZ) THEN
        WIN=0.
        RMT=0.
        WINJ=WIN
      ENDIF
C
      XP=0.D0
      YP=0.D0
      PLABLE(1:NLIM)=.FALSE.
C
      NZAD=5
      NINNE=6
      NIN=20
C
      XMI3D=CH3X0-CH3MX
      XMA3D=CH3X0+CH3MX
      YMI3D=CH3Y0-CH3MY
      YMA3D=CH3Y0+CH3MY
      ZMI3D=CH3Z0-CH3MZ
      ZMA3D=CH3Z0+CH3MZ
      ABSMAX=-1.E30
      ABSMIN=1.D30
      ORDMAX=-1.E30
      ORDMIN=1.D30
C
C  FIND MAXIMA AND MINIMA FOR SCALING OF 3D PLOT WINDOW
C
      LSAVE=PLBOX
      PLBOX=.FALSE.
      CALL PL3D (XMI3D,YMI3D,ZMI3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMI3D,YMI3D,ZMA3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMI3D,YMA3D,ZMI3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMI3D,YMA3D,ZMA3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMA3D,YMI3D,ZMI3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMA3D,YMI3D,ZMA3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMA3D,YMA3D,ZMI3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
      CALL PL3D (XMA3D,YMA3D,ZMA3D,XH,YH)
      ABSMAX=MAX(ABSMAX,XH)
      ABSMIN=MIN(ABSMIN,XH)
      ORDMAX=MAX(ORDMAX,YH)
      ORDMIN=MIN(ORDMIN,YH)
C
C
C
      XNULL=8.
      YNULL=2.
      XCENT=24.
      YCENT=24.
      DELX=ABSMAX-ABSMIN
      DELY=ORDMAX-ORDMIN
C  VERIFY THAT XCENT/YCENT=DELX/DELY
C  IN ORDER TO AVOID RESCALING DURING PLOTTING
C
      YN=XCENT*DELY/DELX
      XN=XCENT
      IF (YN.GT.YCENT) THEN
        XN=YCENT/YN*XCENT
        YN=YCENT
      ENDIF
C
C
      PLBOX=LSAVE
      CALL GRNXTB(1)
      CALL GRSCLC(REAL(XNULL,KIND(1.E0)),REAL(YNULL,KIND(1.E0)),
     .            REAL(XN+XNULL,KIND(1.E0)),
     .            REAL(YN+YNULL,KIND(1.E0)))
      CALL GRSCLV(REAL(ABSMIN,KIND(1.E0)),REAL(ORDMIN,KIND(1.E0)),
     .            REAL(ABSMAX,KIND(1.E0)),REAL(ORDMAX,KIND(1.E0)))
      FAKX=XN/(ABSMAX-ABSMIN)
      FAKY=YN/(ORDMAX-ORDMIN)
C
C  PLOT ADDITIONAL SURFACES
C
C  IF NLTRA, ASSUME THAT THIS SURFACE IS GIVEN IN LOCAL TOROIDAL SYSTEM
C            NO. ILTOR. PLOT IS DONE IN CO-ORDINATE SYSTEM OF CELL ITH
C            WHICH WAS SELECTED BY THE INPUT FLAG ANGLE3
C
      DO 100 I=1,5
        IF (.NOT.PL3A(I)) GOTO 100
        DO 10 IP=1,IPLTA(I)
        DO 10 J=IPLAA(I,IP),IPLEA(I,IP)
          IF (J.GT.NLIMI) GOTO 10
          IF (IGJUM0(J).NE.0) THEN
            IF (TRCPLT) THEN
              WRITE (6,*) 'SURFACE NO. ',J,' OUT'
            ENDIF
            GOTO 10
          ELSE
            IF (NLTRA.AND.ILTOR(J).LE.0) THEN
C  STILL TO BE WRITTEN: BETTER WAY OF IDENTIFYING TOROIDALLY SYMMETRIC
C                       SURFACE
              IF (P3(3,J).GT.1.D50.AND.ZLIMS1(1,J).LT.-1.D15.AND.
     .                                 ZLIMS2(1,J).GT.1.D15) THEN
                LSYMET(J)=.TRUE.
                IF (TRCPLT) THEN
                  WRITE (6,*) 'SURFACE NO. ',J,' TOROIDALLY SYMMETRIC'
                  WRITE (6,*) 'PLOT LATER INTO STANDARD MESH '
                ENDIF
                GOTO 10
              ELSE
                LPERID(J)=.TRUE.
                NJZ=0
                DO 50 IPZ=1,IPLTS(3)
                  IF (IPLAS(3,IPZ).LT.IPLES(3,IPZ)) THEN
                    DO 51 JJ=IPLAS(3,IPZ),IPLES(3,IPZ)-1
                      NJZ=NJZ+1
                      ILT(NJZ)=JJ
51                  CONTINUE
                  ELSE
                    DO 52 JJ=IPLAS(3,IPZ),NTTRA-1
                      NJZ=NJZ+1
                      ILT(NJZ)=JJ
52                  CONTINUE
                    DO 53 JJ=1,IPLES(3,IPZ)-1
                      NJZ=NJZ+1
                      ILT(NJZ)=JJ
53                  CONTINUE
                  ENDIF
50              CONTINUE
                IF (TRCPLT) THEN
                  WRITE (6,*) 'SURFACE NO. ',J,' TOROIDALLY PERODIC'
                  WRITE (6,*) 'PLOT ADD. SURFACE NO. ',J
                  WRITE (6,*) 'INTO ',NJZ, ' TOROIDAL SEGMENTS'
                ENDIF
              ENDIF
            ELSEIF (.NOT.NLTRA.OR.ILTOR(J).GT.0) THEN
              NJZ=1
              ILT(1)=ILTOR(J)
              LSYMET(J)=.FALSE.
              LPERID(J)=.FALSE.
              IF (TRCPLT) THEN
                WRITE (6,*) 'PLOT ADD. SURFACE NO. ',J
              ENDIF
            ENDIF
          ENDIF
C

          DO IJZ=1,NJZ
          IF (NLTRA) WINJ=ZZONE(ILT(IJZ))
C
C** 1 <= RLB < 2 ?
C
          IF (RLB(J).EQ.1..OR.RLB(J).EQ.1.5) THEN
C
C**GEKRUEMMTE FLAECHE ODER  EBENENPAAR ?
            IF (JUMLIM(J).EQ.0) THEN
              CALL FL2O (A0LM(J),A1LM(J),A2LM(J),
     .                   A3LM(J),A4LM(J),A5LM(J),
     .                   A6LM(J),A7LM(J),A8LM(J),
     .                   A9LM(J),MERK,ZX0,ZY0,ZZ0,CX,CY,CZ,
     .                   RZYL,B0,B1,B2,B3,F0,F1,F2,F3,EPS10,NMACH)
              IF (TRCPLT) THEN
                WRITE (6,*) 'FL2O CALLED '
                WRITE (6,*) 'MERK= ',MERK
              ENDIF
C**ZYLINDER: FINDE ACHSE
              IF (MERK.EQ.4) THEN
                IF (TRCPLT) WRITE (6,*) 'CX,CY,CZ,RZYL ',CX,CY,CZ,RZYL
                IF (ABS(CX).GT.1.D-6.AND.
     .              MAX(ABS(CY),ABS(CZ)).LE.1.D-6) THEN
                  T1=(XLIMS1(1,J)-ZX0)/CX
                  T2=(XLIMS2(1,J)-ZX0)/CX
                ELSEIF (ABS(CY).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CZ)).LE.1.D-6) THEN
                  T1=(YLIMS1(1,J)-ZY0)/CY
                  T2=(YLIMS2(1,J)-ZY0)/CY
                ELSEIF (ABS(CZ).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CY)).LE.1.D-6) THEN
                  T1=(ZLIMS1(1,J)-ZZ0)/CZ
                  T2=(ZLIMS2(1,J)-ZZ0)/CZ
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
C  ZYLINDER: GGFLS MEHRERE TEILSTUECKE
                CALL CTQUA (A0LM(J),A1LM(J),A2LM(J),A3LM(J),A4LM(J),
     .                      A5LM(J),A6LM(J),A7LM(J),A8LM(J),A9LM(J),
     .                      XLIMS1(1,J),XLIMS2(1,J),YLIMS1(1,J),
     .                      YLIMS2(1,J),ZLIMS1(1,J),ZLIMS2(1,J),
     .                      RLB(J),RZYL,ZX0,ZY0,ZZ0,PHIAN,PHIEN,IPRT)
                IF (TRCPLT) THEN
                  WRITE (6,*) ' T1,T2 ',T1,T2
                  WRITE (6,*) ' IPRT ',IPRT
                  WRITE (6,*) ' PHIAN,PHIEN ',(PHIAN(JP),PHIEN(JP),
     .                                         JP=1,IPRT)
                ENDIF
                DO 109 IPR=1,IPRT
                  PHA=PHIAN(IPR)
                  PHE=PHIEN(IPR)
                  CALL ZYLIND (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                         RZYL,NZAD,NINNE,NIN,ILCOL(J),
     .                         IGFIL(J).NE.0,
     .                         J,0,AL,0,AR,PHA,PHE)
109             CONTINUE
C**KEGEL: BISLANG NUR EIN STUECK MOEGLICH. FINDE ACHSE
              ELSEIF (MERK.EQ.8) THEN
                IF (TRCPLT) WRITE (6,*) 'CX,CY,CZ,RZYL ',CX,CY,CZ,RZYL
                IF (ABS(CX).GT.1.D-6.AND.
     .              MAX(ABS(CY),ABS(CZ)).LE.1.D-6) THEN
                  T1=(XLIMS1(1,J)-ZX0)/CX
                  T2=(XLIMS2(1,J)-ZX0)/CX
                ELSEIF (ABS(CY).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CZ)).LE.1.D-6) THEN
                  T1=(YLIMS1(1,J)-ZY0)/CY
                  T2=(YLIMS2(1,J)-ZY0)/CY
                ELSEIF (ABS(CZ).GT.1.D-6.AND.
     .                  MAX(ABS(CX),ABS(CY)).LE.1.D-6) THEN
                  T1=(ZLIMS1(1,J)-ZZ0)/CZ
                  T2=(ZLIMS2(1,J)-ZZ0)/CZ
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
                CALL CONE (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                     RZYL,NZAD,NINNE,NIN,ILCOL(J),
     .                     IGFIL(J).NE.0,J,0,AL,0,AR)
C**KUGEL,ELLIPSOID: BISLANG NUR EIN STUECK MOEGLICH
C             ELSEIF (MERK.EQ.8) CALL SPHERE (ZX0,ZY0,ZZ0,HX,HY,HZ,
C    .                                   NZAD,NINNE,NIN,ILCOL(J),
c    .                                   IGFIL(J).NE.0,J,0,AL,0,AR)
C**KUGEL, ELLIPSOID
              ELSEIF (MERK == 13) THEN
                CALL ELLIPSOID (ZX0,ZY0,ZZ0,CX,CY,CZ,XLIMS1(1,J),
     .               YLIMS1(1,J),ZLIMS1(1,J),XLIMS2(1,J),YLIMS2(1,J),
     .               ZLIMS2(1,J),RLB(J),ILCOL(J),5,5,5)
C**PAAR VON EBENEN (ODER EINE DOPPELEBENE)
              ELSEIF (MERK.EQ.1.OR.MERK.EQ.2.OR.MERK.EQ.3) THEN
                IF (TRCPLT) WRITE (6,*) 'B0,B1,B2,B3 ',B0,B1,B2,B3
                CALL PLANE (B0,B1,B2,B3,RLB(J),9,EPS10,
     .                      ALIMS,XLIMS,YLIMS,ZLIMS,
     .                      ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                             XLIMS2,YLIMS2,ZLIMS2,
     .                             XLIMS3,YLIMS3,ZLIMS3,
     .                      ILCOL(J),IGFIL(J).NE.0,J)
                IF (MERK.NE.1) THEN
                  IF (TRCPLT) WRITE (6,*) 'F0,F1,F2,F3 ',F0,F1,F2,F3
                  CALL PLANE (F0,F1,F2,F3,RLB(J),9,EPS10,
     .                        ALIMS,XLIMS,YLIMS,ZLIMS,
     .                        ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                               XLIMS2,YLIMS2,ZLIMS2,
     .                               XLIMS3,YLIMS3,ZLIMS3,
     .                        ILCOL(J),IGFIL(J).NE.0,J)
                ENDIF
              ELSE
                PLABLE(J)=.TRUE.
                GOTO 10
              ENDIF
C**EINE EBENE
            ELSEIF (JUMLIM(J).NE.0) THEN
              CALL PLANE (A0LM(J),A1LM(J),A2LM(J),A3LM(J),RLB(J),9,
     .              EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                    ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                           XLIMS2,YLIMS2,ZLIMS2,
     .                           XLIMS3,YLIMS3,ZLIMS3,
     .                    ILCOL(J),IGFIL(J).NE.0,J)
            ENDIF
C
C**RLB >= 3 ? EIN EBENENSTUECK, DURCH POLYGON BEGRENZT
C
          ELSEIF (RLB(J).GT.2.) THEN
C
            CALL PRLLO(P1(1,J),P2(1,J),P3(1,J),P4(1,J),P5(1,J),
     .                 ILCOL(J),IGFIL(J).NE.0)
C
C**RLB < 0 ? ERST EINIGE OPTIONEN VORHANDEN, REST: CALL PLTUSR
C
          ELSEIF (RLB(J).LT.0) THEN
C
            IF (JUMLIM(J).EQ.0) THEN
              CALL FL2O (A0LM(J),A1LM(J),A2LM(J),
     .                   A3LM(J),A4LM(J),A5LM(J),
     .                   A6LM(J),A7LM(J),A8LM(J),
     .                   A9LM(J),MERK,ZX0,ZY0,ZZ0,CX,CY,CZ,
     .                   RZYL,B0,B1,B2,B3,F0,F1,F2,F3,EPS10,NMACH)
              IF (TRCPLT) THEN
                WRITE (6,*) 'FL2O CALLED '
                WRITE (6,*) 'MERK= ',MERK
              ENDIF
C
C**PAAR VON EBENEN ODER DOPPELEBENE ?
              IF (MERK.LE.3) THEN
                IF (ISCN(J).EQ.0) THEN
                  CALL PLANE (B0,B1,B2,B3,RLB(J),
     .                        9,EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                        ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                               XLIMS2,YLIMS2,ZLIMS2,
     .                               XLIMS3,YLIMS3,ZLIMS3,
     .                        ILCOL(J),IGFIL(J).NE.0,J)
                  CALL PLANE (F0,F1,F2,F3,RLB(J),
     .                        9,EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                        ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                               XLIMS2,YLIMS2,ZLIMS2,
     .                               XLIMS3,YLIMS3,ZLIMS3,
     .                        ILCOL(J),IGFIL(J).NE.0,J)
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
C**ZYLINDER ?
              ELSEIF (MERK.EQ.4) THEN
C**ZYLINDER BEGRENZT DURCH MAXIMAL 9 EBENEN
               IF (ISCN(J).EQ.0) THEN
                 CALL ZYLPLN (ZX0,ZY0,ZZ0,CX,CY,CZ,RZYL,J,NZAD,NINNE,
     .                        NIN)
C**ZYLINDER BEGRENZT DURCH MAXIMAL EINE FLAECHE ZWEITER ORDNUNG
               ELSEIF (ILIN(J).EQ.0.AND.ISCN(J).EQ.1) THEN
                 IB=1
                 CALL FL2O (ALIMS0(IB,J),XLIMS1(IB,J),YLIMS1(IB,J),
     .                      ZLIMS1(IB,J),XLIMS2(IB,J),YLIMS2(IB,J),
     .                      ZLIMS2(IB,J),XLIMS3(IB,J),YLIMS3(IB,J),
     .                      ZLIMS3(IB,J),MERK2,ZX0B,ZY0B,ZZ0B,CXB,CYB,
     .                      CZB,RZYLB,B0B,B1B,B2B,B3B,F0B,F1B,F2B,F3B,
     .                      EPS10,NMACH)
                 IF (TRCPLT) THEN
                   WRITE (6,*) 'FL2O CALLED '
                   WRITE (6,*) 'MERK2= ',MERK2
                 ENDIF
C**ZYLINDER BEGRENZT VON 2 EBENEN
                 IF (MERK2.LE.3) THEN
                   AL(1)=B0B
                   AL(2)=B1B
                   AL(3)=B2B
                   AL(4)=B3B
                   AR(1)=F0B
                   AR(2)=F1B
                   AR(3)=F2B
                   AR(4)=F3B
                   CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AL,4,TA,TD,LERR1)
                   CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AR,4,TB,TD,LERR2)
                   IF (LERR1.OR.LERR2) THEN
                     IF (TRCPLT)
     .               WRITE (6,*)' FEHLER IN BERANDUNG VON FLAECHE ',J
                     PLABLE(J)=.TRUE.
                     GOTO 10
                   ENDIF
                   IF (TA.LT.TB) THEN
                     T1=TA-2.*RZYL
                     T2=TB+2.*RZYL
                   ELSE
                     T1=TB-2.*RZYL
                     T2=TA+2.*RZYL
                     DO 15 II=1,4
                       TD=AL(II)
                       AL(II)=AR(II)
                       AR(II)=TD
15                   CONTINUE
                   ENDIF
                   CALL ZYLIND (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                          RZYL,NZAD,NINNE,NIN,
     .                          ILCOL(J),IGFIL(J).NE.0,J,4,AL,4,AR,
     .                          0._DP,360._DP)
C**ZYLINDER BEGRENZT VON ECHT GEKRUEMMTEN FLAECHE 2TER ORDNUNG
                ELSEIF (MERK2.GE.4) THEN
                  IB=1
                  AL(1)=ALIMS0(IB,J)
                  AL(2)=XLIMS1(IB,J)
                  AL(3)=YLIMS1(IB,J)
                  AL(4)=ZLIMS1(IB,J)
                  AL(5)=XLIMS2(IB,J)
                  AL(6)=YLIMS2(IB,J)
                  AL(7)=ZLIMS2(IB,J)
                  AL(8)=XLIMS3(IB,J)
                  AL(9)=YLIMS3(IB,J)
                  AL(10)=ZLIMS3(IB,J)
                  CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AL,10,TA,TB,LERR1)
                  IF (LERR1) THEN
                    WRITE (6,*)' FEHLER IN DER BERANDUNG VON FLAECHE',J
                    PLABLE(J)=.TRUE.
                    GOTO 10
                  ENDIF
                  IF (TA.LT.TB) THEN
                    T1=TA-2.*RZYL
                    T2=TB+2.*RZYL
                  ELSE
                    T1=TB-2.*RZYL
                    T2=TA+2.*RZYL
                  ENDIF
                  DO 18 K=1,10
                    AR(K)=AL(K)
18                CONTINUE
                  CALL ZYLIND (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                      RZYL,NZAD,NINNE,NIN,
     .                      ILCOL(J),IGFIL(J).NE.0,
     .                      J,10,AL,10,AR,0._DP,360._DP)
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
               ENDIF
C**KUGEL, ELLIPSOID
              ELSEIF (MERK == 13) THEN
                CALL ELLIPSOID (ZX0,ZY0,ZZ0,CX,CY,CZ,XLIMS1(1,J),
     .               YLIMS1(1,J),ZLIMS1(1,J),XLIMS2(1,J),YLIMS2(1,J),
     .               ZLIMS2(1,J),RLB(J),ILCOL(J),5,5,5)
              ELSE
                PLABLE(J)=.TRUE.
                GOTO 10
              ENDIF
C
C**EBENE MIT RLB.LT.0 OPTION
C
            ELSEIF (JUMLIM(J).NE.0) THEN
C
C**EBENE BEGRENZT DURCH ANDERE EBENEN
              IF (ISCN(J).EQ.0) THEN
                CALL PLANE (A0LM(J),A1LM(J),A2LM(J),A3LM(J),RLB(J),
     .                      9,EPS10,ALIMS,XLIMS,YLIMS,ZLIMS,
     .                      ALIMS0,XLIMS1,YLIMS1,ZLIMS1,
     .                             XLIMS2,YLIMS2,ZLIMS2,
     .                             XLIMS3,YLIMS3,ZLIMS3,
     .                      ILCOL(J),IGFIL(J).NE.0,J)
              ELSEIF (ILIN(J).EQ.0) THEN
C**EBENE BEGRENZT DURCH EINEN ODER MEHRERE ZYLINDER?
                IB=0
20              IB=IB+1
                CALL FL2O (ALIMS0(IB,J),XLIMS1(IB,J),YLIMS1(IB,J),
     .                     ZLIMS1(IB,J),XLIMS2(IB,J),YLIMS2(IB,J),
     .                     ZLIMS2(IB,J),XLIMS3(IB,J),YLIMS3(IB,J),
     .                     ZLIMS3(IB,J),MERK2,ZX0,ZY0,ZZ0,CX,CY,CZ,
     .                     RZYL,B0,B1,B2,B3,F0,F1,F2,F3,EPS10,NMACH)
                IF (TRCPLT) THEN
                  WRITE (6,*) 'FL2O CALLED '
                  WRITE (6,*) 'MERK2= ',MERK2
                ENDIF
                IF (MERK2.EQ.4) THEN
                  AL(1)=A0LM(J)
                  AL(2)=A1LM(J)
                  AL(3)=A2LM(J)
                  AL(4)=A3LM(J)
                  CALL SECQUA (ZX0,ZY0,ZZ0,CX,CY,CZ,AL,4,TA,TD,LERR1)
                  IF (LERR1) THEN
                    WRITE (6,*)' FEHLER IN DER BERANDUNG VON FLAECHE',J
                    PLABLE(J)=.TRUE.
                    GOTO 10
                  ENDIF
                  T1=TA-2.*RZYL
                  T2=TA+4.*RZYL
                  NP=1
                  CALL ZYLIND (ZX0,ZY0,ZZ0,CX,CY,CZ,T1,T2,
     .                         RZYL,NP,NINNE,NIN,
     .                         ILCOL(J),IGFIL(J).NE.0,J,4,AL,0,AR,
     .                         0._DP,360._DP)
                ELSE
                  PLABLE(J)=.TRUE.
                  GOTO 10
                ENDIF
                IF (IB.LT.ISCN(J)) GOTO 20
C**EBENE BEGRENZT DURCH ALLE ANDERE OPTIONEN
              ELSE
                PLABLE(J)=.TRUE.
                GOTO 10
              ENDIF
C
            ENDIF
          ENDIF
C
C END NJZ LOOP
          ENDDO
C
10      CONTINUE
100   CONTINUE
C
C
C  PLOTTE DIEJENIGEN FLAECHEN, DIE NICHT AUTOMATISCH
C  MOEGLICH WAREN, IN DER USER-SUPPLIED ROUTINE PLTUSR
      I1=0
      DO 200 J=1,NLIMI
        IF (PLABLE(J)) THEN
          IF (NLTRA) WINJ=ZZONE(ILTOR(J))
          CALL PLTUSR(PLABLE(J),J)
        ENDIF
C
        IF (PLABLE(J)) THEN
          IF (I1.EQ.0) THEN
            WRITE (6,*) 'MESSAGE FROM SUBR. PLT3D :'
            WRITE (6,*) 'SURFACE NO. J COULD NOT BE PLOTTED'
            I1=1
          ENDIF
          WRITE (6,*) 'J= ',J
        ENDIF
200   CONTINUE
C
C  PLOT SURFACES OF STANDARD MESH
C
500   IF (.NOT.(PL3S(1).OR.PL3S(2).OR.PL3S(3))) GOTO 10000
C
C  TOROIDAL GRID
C
      DO 3000 IPZ=1,IPLTS(3)
C
        NJZ=0
        IF (IPLAS(3,IPZ).LT.IPLES(3,IPZ)) THEN
          DO 3001 J=IPLAS(3,IPZ),IPLES(3,IPZ)
            NJZ=NJZ+1
            IF (J.EQ.1.OR.(NLTOR.OR.NLTRA)) THEN
              ZPLOT(NJZ)=ZSURF(J)
            ELSEIF (J.EQ.2) THEN
              ZPLOT(NJZ)=ZAA
            ENDIF
3001      CONTINUE
        ELSE
          DO 3002 J=IPLAS(3,IPZ),NTTRA
            NJZ=NJZ+1
            ZPLOT(NJZ)=ZSURF(J)
3002      CONTINUE
          DO 3003 J=2,IPLES(3,IPZ)
            NJZ=NJZ+1
            ZPLOT(NJZ)=ZSURF(J)
3003      CONTINUE
        ENDIF
C
        DO 3100 IZ=1,NJZ
          CALL GRNWPN(2)
          PHI=ZPLOT(IZ)
C
C  PHI = CONST , PLOT POLOIDAL CROSS SECTION AT TOROIDAL POSITION PHI
C
          IST=5
          IS=0
C
          DO 1000 IBR=1,IPLTS(1)
C
            IS=0
            IST=10
C
            DO 1100 IR=IPLAS(1,IBR),IPLES(1,IBR)
              IF (TRCPLT) WRITE (6,*) 'PLOT RAD. STAND. SURFACE NO. ',IR
C
C NR: POINTS TO BE PLOTTED ON RADIAL SURFACE IR
C
              IF (NLCRC.OR.NLELL.OR.NLTRI) THEN
                DM=0.
                RS=RSURF(IR)
                EP=EP1(IR)
                EL=ELL(IR)
                TR=TRI(IR)
                CALL PLGELR(RS,EP,EL,TR,DM,100,XX,YY,NR,PSURF,NP2ND)
C
                DO 1120 J=1,NR
                  Y=YY(J)
                  IF (NLTRA) THEN
                    RR=XX(J)+RMTOR
                    X=RR*COS(PHI)
                    Z=RR*SIN(PHI)
                    CALL TORLOC(WIN,RMT,X,Z)
                    WINJ=WIN
                  ELSEIF (NLTRZ) THEN
                    X=XX(J)
                    Z=PHI
                  ENDIF
C
                  CALL PL3D(X,Y,Z,XP(J),YP(J))
1120            CONTINUE
                do 1130 jj=1,nr
                  xps(jj)=xp(jj)
                  yps(jj)=yp(jj)
1130            continue
                CALL GRLN(XPS,YPS,NR)
C
C
              ELSEIF (NLPLG) THEN
C
                NR=0
                DO 1160 K=1,NPPLG
                  IAN=NPOINT(1,K)
                  IEN=NPOINT(2,K)
                  KIN=NR+1
                  DO 1165 J=IAN,IEN
                    IF (NR.GE.NPLY) THEN
                      WRITE (6,*) 'FROM PLT3D: NOT ENOUGH STORAGE   '
                      WRITE (6,*) 'INCREASE PARAMETER NPLY '
                      GOTO 1165
                    ENDIF
                    NR=NR+1
                    X=XPOL(IR,J)
                    Y=YPOL(IR,J)
                    IF (NLTRA) THEN
                      RR=X+RMTOR
                      X=RR*COS(PHI)
                      Z=RR*SIN(PHI)
                      CALL TORLOC(WIN,RMT,X,Z)
                      WINJ=WIN
                    ELSEIF (NLTRZ) THEN
                      Z=PHI
                    ENDIF
                    CALL PL3D(X,Y,Z,XP(NR),YP(NR))
1165              CONTINUE
                  do 1162 jj=1,nr
                    xps(jj)=xp(jj)
                    yps(jj)=yp(jj)
1162              continue
                  CALL GRLN (XPS,YPS,NR)
1160            CONTINUE
C
C
              ELSE
C  TO BE WRITTEN
              ENDIF
C
C  RADIAL SURFACE IR PLOTTED, NR POINTS
C  NEXT: SAVE CO-ORDINATES
C
              IF (IS+NR/IST+1.LE.NPLY) THEN
                DO 1200 J=1,NR,IST
                  IS=IS+1
                  XSAVE(IS,IZ)=XP(J)
                  YSAVE(IS,IZ)=YP(J)
1200            CONTINUE
C  LAST POINT:
                IF (MOD(NR,IST).NE.0) THEN
                  IS=IS+1
                  XSAVE(IS,IZ)=XP(NR)
                  YSAVE(IS,IZ)=YP(NR)
                ENDIF
              ELSE
                WRITE(6,*) ' STORAGE EXCEEDED IN PLT3D, IR= ',IR
              ENDIF
C
1100        CONTINUE
C
C  ALL RADIAL SURFACES IN BLOCK IBR DONE
C
1000      CONTINUE
C
C  ALL BLOCKS FOR RADIAL SURFACES DONE
C
C  ARE THERE TOROIDALLY (OR Z) SYMMETRIC ADDITIONAL SURFACES
C
          ISSTD=IS
          CALL GRNWPN(1)
          LZR=.FALSE.
C
          CH2X0S=CH2X0
          CH2Y0S=CH2Y0
          CH2MXS=CH2MX
          CH2MYS=CH2MY
          CH2X0=CH3X0
          CH2Y0=CH3Y0
          CH2MX=CH3MX
          CH2MY=CH3MY
          DO 1300 IBA=1,5
            IF (.NOT.PL3A(IBA)) GOTO 1300
            DO 1310 IA=1,IPLTA(IBA)
              DO 1320 J=IPLAA(IBA,IA),IPLEA(IBA,IA)
                IF (IGJUM0(J).NE.0) GOTO 1320
                IF (.NOT.LSYMET(J).OR.J.GT.NLIMI) GOTO 1320
                PLT1=PLCUT(1)
                PLT2=PLCUT(2)
                PLT3=PLCUT(3)
                PLCUT(1)=.FALSE.
                PLCUT(2)=.FALSE.
                PLCUT(3)=.TRUE.
                CALL PLTADD(J,J)
                PLCUT(1)=PLT1
                PLCUT(2)=PLT2
                PLCUT(3)=PLT3
                NA=0
C  AT PRESENT: ONLY FIRST AND LAST POINT, IF STRAIGHT LINE
C              OR 10 POINTS, IF CURVED LINE
                ISTP=MAX(1,(INSTOR-1)/10)
                IF (JUMLIM(J).GT.0) ISTP=INSTOR-1
                IF (INSTOR.LT.2) GOTO 1320
                CUR => FIRST_POINT
                DO WHILE(ASSOCIATED(CUR))
                  NA=NA+1
                  X=CUR%XPL2D
                  Y=CUR%YPL2D
                  IF (NLTRA) THEN
                    RR=X+RMTOR
                    X=RR*COS(PHI)
                    Z=RR*SIN(PHI)
                    CALL TORLOC(WIN,RMT,X,Z)
                    WINJ=WIN
                  ELSEIF (NLTRZ) THEN
                    Z=PHI
                  ENDIF
                  CALL PL3D(X,Y,Z,XP(NA),YP(NA))
                  IF (CUR%NPL2D.EQ.0) CALL GRJMP (
     .               REAL(XP(NA),KIND(1.E0)),
     .               REAL(YP(NA),KIND(1.E0)))
                  IF (CUR%NPL2D.EQ.1) CALL GRDRW (
     .               REAL(XP(NA),KIND(1.E0)),
     .               REAL(YP(NA),KIND(1.E0)))
C
                  IS=IS+1
                  XSAVE(IS,IZ)=XP(NA)
                  YSAVE(IS,IZ)=YP(NA)
                  DO ID=1,ISTP
                    IF (ASSOCIATED(CUR)) CUR => CUR%NXTPNT
                  END DO
                END DO
1320          CONTINUE
1310        CONTINUE
1300      CONTINUE
          CH2X0=CH2X0S
          CH2Y0=CH2Y0S
          CH2MX=CH2MXS
          CH2MY=CH2MYS
C
3100    CONTINUE
C
C   LINES OF CONSTANT POLOIDAL POSITION
C
        CALL GRDSH(0.2,0.5,0.2)
        CALL GRNWPN(2)
        DO 2000 J=1,IS
          IF (J.GT.ISSTD) CALL GRNWPN(1)
          CALL GRJMP (REAL(XSAVE(J,1),KIND(1.E0)),
     .                REAL(YSAVE(J,1),KIND(1.E0)))
          DO 2000 IZ=2,NJZ
            CALL GRDRW(REAL(XSAVE(J,IZ),KIND(1.E0)),
     .                 REAL(YSAVE(J,IZ),KIND(1.E0)))
2000    CONTINUE
        CALL GRDSH(1.,0.,1.)
        CALL GRNWPN(1)
C
3000  CONTINUE
C
10000 CONTINUE
C
C  BESCHRIFTUNG
C
      XH=(ABSMIN+ABSMAX)/2.
      YH=ORDMAX+2./FAKY
      CALL GRTXT (REAL(XH,KIND(1.E0)),REAL(YH,KIND(1.E0)),27,
     .            'CHECK OF GEOMETRICAL INPUT:')
      YH=YH-0.5/FAKY
      CALL GRTXT (REAL(XH,KIND(1.E0)),REAL(YH,KIND(1.E0)),72,TXTRUN)
      YH=YH-0.75/FAKY
      DO 11000 J=1,5
        IF (PL3A(J)) CALL GRTXT(REAL(XH,KIND(1.E0)),REAL(YH,KIND(1.E0)),
     .                          16,TEXTLA(J))
        IF (PL3A(J)) YH=YH-0.5/FAKY
11000 CONTINUE
C
C  RETURN CO-ORDINATES FOR PLOTS OF PARTICLE TRACKS
      XR=ABSMIN-6./FAKX
      YR=ORDMAX-0./FAKY
      RETURN
      END
