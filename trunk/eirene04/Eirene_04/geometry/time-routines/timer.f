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
     .           MPOL, IPOLGO, LEARC2, IP, ICELLR, MSAVE
      LOGICAL :: LCUT(N2NDPLG)
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
C
C     IF (NLTRC) THEN
C       CALL LEER(1)
C       WRITE (6,*) 'TIMER, INIT.: NRCELL,NLSRFX,MRSURF,TT,TL'
C       WRITE (6,*)                NRCELL,NLSRFX,MRSURF,TT,TL
C       WRITE (6,*) '              NPCELL,NLSRFY,MPSURF'
C       WRITE (6,*)                NPCELL,NLSRFY,MPSURF
C     ENDIF
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
C     IF (NLTRC) THEN
C       WRITE (6,*) 'TIMER, OUT: PT,MRSURF,NINCX'
C       WRITE (6,*)              PT,MRSURF,NINCX
C     ENDIF
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
C     IF (NLTRC) THEN
C       WRITE (6,*) 'NJUMP,NINCX,NRCELL '
C       WRITE (6,*)  NJUMP,NINCX,NRCELL
C     ENDIF
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
        WRITE (6,*) 'NRCELL.LE.0 IN TIMER: NOT READY '
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
C     IF (NLTRC) WRITE (6,*) 'IPOLGN,NRCELL FROM TIMER ',IPOLGN,NRCELL
      IF (NRI.EQ.0) THEN
        RETURN
      ELSE
        WRITE (6,*) ' OPTION NRI > 0 NOT READY in TIMER '
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
C       IF (NLTRC) WRITE (6,*) ' TIMER: IN NEIGHBOR PART'
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
C       IF (NLTRC) WRITE (6,*) ' IR,IP,LNGB1,LNGB2,LNGB3,LNGB4',
C    .                           IR,IP,LNGB1,LNGB2,LNGB3,LNGB4
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
C       IF (NLTRC) WRITE (6,*) ' T1,T2,T3,T4 ',T1,T2,T3,T4
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
C       IF (NLTRC) WRITE (6,*) ' LCT1,LCT2,LCT3,LCT4 ',
C    .                           LCT1,LCT2,LCT3,LCT4
C       IF (NLTRC) WRITE (6,*) ' PT1,PT2,PT3,PT4 ',PT1,PT2,PT3,PT4
        T=MAX(PT1,PT2,PT3,PT4)
C       IF (NLTRC) WRITE (6,*) ' T = ',T
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
            WRITE (6,*) ' NO INTERSECTION IN TIMER. TRY ONCE AGAIN '
            MMSURF=MSURF
            IF (MSURF.GT.NLIM) MMSURF=-MSURF+NLIM
            WRITE (6,*) 'NPANU, MSURF ',NPANU,MMSURF
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
          WRITE (6,*) ' ERROR: NO INTERSECTION FOUND IN TIMER '
          WRITE (6,*) ' NPANU: ',NPANU
          MRSURF=0
          PT=1.D30
          NINCX=0
          ICELLR=0
          IPOLGN=IP
        ENDIF
        NJUMP=1
C       IF (NLTRC) THEN
C         WRITE (6,*) ' PT,MRSURF,NINCX,IPOLGN ',
C    .                  PT,MRSURF,NINCX,IPOLGN
C         WRITE (6,*) 'NCOUP ',NCOUP
C         DO ICOUP=1,NCOUP
C           WRITE (6,*) 'ICOUP,ALPD(ICOUP) ',ICOUP,ALPD(ICOUP)
C         ENDDO
C       ENDIF
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
C         IF (NLTRC) THEN
C           WRITE (6,*) ' SURFACE NO. = ',I
C           WRITE (6,*) ' TIMPOL ',(TIMPOL(I,J),J=1,NTIM(I))
C           WRITE (6,*) ' IIMPOL ',(IIMPOL(I,J),J=1,NTIM(I))
C         ENDIF
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
C     IF (NLTRC) WRITE (6,*) 'PT,MRSURF,IPOLGN,NINCX ',
C    .                        PT,MRSURF,IPOLGN,NINCX
      RETURN
C
C
8000  CONTINUE
C
      IF (LEVGEO.GT.4) GOTO 10000
C
C   FINITE ELEMENT DISCRETISATION
C
C     IF (NLTRC) THEN
C       WRITE (6,*) ' TIMER,NJUMP,NRCELL ',NJUMP,NRCELL
C       WRITE (6,*) '       NLSRFX,IPOLG ',NLSRFX,IPOLG
C     ENDIF
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
C       IF (NLTRC) THEN
C         WRITE (6,*) ' IZELL,IPOLGO ',IZELL,IPOLGO
C         WRITE (6,*) ' X0,Y0 ',X0,Y0
C         WRITE (6,*) ' VELX,VELY ',VELX,VELY
C         CALL LEER(1)
C       ENDIF
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
C  NUMMER DER GESCHNITTENEN SEITE DES ALTEN DREIECKS
                NTIM(IZELL) = J
C               IF (NLTRC) THEN
C                 WRITE (6,*) 'IZELL,XX,YY,J ',IZELL,XX,YY,J
C               ENDIF
C  ZELLENNUMMER DES NEUEN DREIECKS
                IZELLO=IZELL
                IZELL = NCHBAR(J,IZELLO)
                IF (IZELL .EQ. 0) GOTO 8050
C  SEITENNUMMER DES NEUEN DREIECKS
                IPOLGO = NSEITE(J,IZELLO)
C               IF (NLTRC) THEN
C                 WRITE(6,*) 'GEHE IN ZELLE ',IZELL,' TM= ',TM
C                 WRITE(6,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
C               ENDIF
                GOTO 8050
              ENDIF
            ENDIF
          ENDIF
8010    CONTINUE
C
C       IF (NLTRC) WRITE (6,*) ' NO INTERSECTION FOUND IN TRIANGLE '
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
C     IF (NLTRC) WRITE (6,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
C    .                         NRCELL,MRSURF,PT,NINCX,IPOLGN

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
C       IF (NLTRC) THEN
C         WRITE (6,*) ' IZELL,IPOLGO ',IZELL,IPOLGO
C         WRITE (6,*) ' X0,Y0,Z0 ',X0,Y0,Z0
C         WRITE (6,*) ' VELX,VELY,VELZ ',VELX,VELY,VELZ
C       ENDIF

        XX = X0
        YY = Y0
        ZZ = Z0
        TM=0.
        IF (IZELL.EQ.0) THEN
          WRITE (6,*) ' IZELL = 0 IN TIMER '
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
            PNORMI=RINCRC(J,IZELL)
            A(1:3,1) = (/ VTETX(I1,IZELL), VTETY(I1,IZELL),
     .                    VTETZ(I1,IZELL)/) * SIG1 * PNORMI
            A(1:3,2) = (/ VTETX(I2,IZELL), VTETY(I2,IZELL),
     .                    VTETZ(I2,IZELL)/) * SIG2 * PNORMI
            A(1:3,3) = (/ -VELX, -VELY, -VELZ /) * PNORMI
            B(1:3) = (/ XX-XTETRA(NTECK(ITSIDE(1,J),IZELL)),
     .                  YY-YTETRA(NTECK(ITSIDE(1,J),IZELL)),
     .                  ZZ-ZTETRA(NTECK(ITSIDE(1,J),IZELL)) /) * PNORMI
cpb            DET = SARRUS(A)
            DET = A(1,1) * A(2,2) * A(3,3) 
     .          + A(1,2) * A(2,3) * A(3,1)
     .          + A(1,3) * A(2,1) * A(3,2)
     .          - A(3,1) * A(2,2) * A(1,3)
     .          - A(3,2) * A(2,3) * A(1,1)
     .          - A(3,3) * A(2,1) * A(1,2)
            IF (ABS(DET) < 1.D-5) CYCLE
            AB = A
            AB(:,1) = B
cpb            XMU = SARRUS(AB)/DET
            XMU = ( AB(1,1) * AB(2,2) * AB(3,3) 
     .            + AB(1,2) * AB(2,3) * AB(3,1)
     .            + AB(1,3) * AB(2,1) * AB(3,2)
     .            - AB(3,1) * AB(2,2) * AB(1,3)
     .            - AB(3,2) * AB(2,3) * AB(1,1)
     .            - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
            IF (XMU .GE.0.D0 .AND. XMU .LE.1.D0) THEN
              AB = A
              AB(:,2) = B
cpb              XETA = SARRUS(AB)/DET
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
cpb                T = SARRUS(AB)/DET
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
                  TM = TM + T
                  TIMTET = TM
C  NUMMER DER GESCHNITTENEN SEITE DES ALTEN DREIECKS
                  NTIMT = J
C                 IF (NLTRC) THEN
C                   WRITE (6,*) 'IZELL,XX,YY,ZZ,J ',IZELL,XX,YY,ZZ,J
C                 ENDIF
C  ZELLENNUMMER DES NEUEN DREIECKS
                  IZELLO=IZELL
                  IPOLGOO=IPOLGO
                  IZELL = NTBAR(J,IZELLO)
                  PT = TM
                  IF (IZELL .EQ. 0) GOTO 11050
C  SEITENNUMMER DES NEUEN DREIECKS
                  IPOLGO = NTSEITE(J,IZELLO)
C                 IF (NLTRC) THEN
C                   WRITE(6,*) 'GEHE IN ZELLE ',IZELL,' TM= ',TM
C                   WRITE(6,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
C                 ENDIF
                  GOTO 11050
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        END DO
C
C       IF (NLTRC) WRITE (6,*) ' NO INTERSECTION FOUND IN TETRAHEDRON '
        PT=1.D30
        NINCX=0
C  DURCHSUCHE NACHBARDREIECK
        IF (NLSRFX) THEN
          IZELL=NTBAR(IPOLG,NRCELL)
          IPOLGO=NTSEITE(IPOLG,NRCELL)
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
C       IF (NLTRC) WRITE (6,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
C    .                           NRCELL,MRSURF,PT,NINCX,IPOLGN

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
C         IF (NLTRC) THEN
C           WRITE(6,*) 'OUTSIDE; GEHE IN ZELLE ',NTMZ,' TM= ',TM
C           WRITE(6,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
C         ENDIF
        END IF

        NLSRFX=.FALSE.
        NJUMP=1
        MRSURF=NTMZ
        IPOLGN=NTMS
        NINCX=NTMZ-NRCELL
C       IF (NLTRC) WRITE (6,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
C    .                           NRCELL,MRSURF,PT,NINCX,IPOLGN
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
            WRITE (6,*) ' ERROR IN TIMER, NRMSRF WRONG '
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
        WRITE (6,*) 'NEWCEL=0, EXIT FROM MESH '
        CALL EXIT_OWN(1)
      ENDIF
      NLSRFX=.FALSE.
      RETURN

C
9999  CONTINUE
      WRITE (6,*) 'ERROR IN TIMER, EXIT CALLED AT 9999'
      WRITE (6,*) 'NPANU ', NPANU
      CALL EXIT_OWN(1)
      END
