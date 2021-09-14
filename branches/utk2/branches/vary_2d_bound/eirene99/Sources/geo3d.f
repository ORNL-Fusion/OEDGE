C
C  FULL EIRENE GEOMETRY BLOCK  (GEO3D)
C
C
      SUBROUTINE TIMER (PT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CLOGAU'
      INCLUDE 'CPOLYG'
      INCLUDE 'CCONA'
      INCLUDE 'CLGIN'
      INCLUDE 'CUPD'
      INCLUDE 'CTRIG'
c slmod begin - grid - not tr
      INTEGER          CheckCell,PolPos,CalcInter
      LOGICAL          funnyoutput
      INTEGER          AP,VP1,ILOOP
      DATA funnyoutput /.TRUE./
c slmod end
C
      LOGICAL LCUT(N2NDPLG),
     L        LNGB1,LNGB2,LNGB3,LNGB4,LCT1,LCT2,LCT3,LCT4
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
c slmod begin - tr
      ILOOP=0
c slmod end
      IRCELL=NRCELL
      PT=1.D30
      IF (ABS(VELZ).EQ.1.D0) THEN
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
      IF (NRCELL.GT.0) THEN
        MRSURF=NRCELL
        IF (NINCX.EQ.1) MRSURF=MRSURF+1
        DXA=RSURF(MRSURF)-XA
        PT=DXA/(VELX+EPS60)
      ELSE
C  TRY TO FIND REENTRY SURFACE. CHECK ONLY NON DEFAULT RADIAL SURFACES
        MRSURF=0
        DO 10 IR=1,NR1ST
          IF (INMP1I(IR,0,0).NE.0) THEN
            DXA=RSURF(IR)-XA
            TST=DXA/(VELX+EPS60)
            IF (TST.LE.0..OR.TST.GT.PT) GOTO 10
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
      IF (NJUMP.EQ.0) THEN
        IF (NLSRFX) THEN
          ZC1=RQ(MRSURF)
        ELSE
          ZC1=X0*X0+Y0*Y0
        ENDIF
        ZA=VELX*VELX+VELY*VELY
        ZB=X0*VELX+Y0*VELY
        ZB2=ZB*ZB
        ZAB=-ZB/(ZA+EPS60)
        ZAB2=ZA/(ZB2+EPS60)
C
C  TEST FOR DIRECTION
        NINCX=1
        IF (ZB.LT.0) NINCX=-NINCX
      ENDIF
C
      IF (NRCELL.GT.0) THEN
C  PARTICLE INSIDE STANDARD MESH
        NRI=0
        MRSURF=NRCELL
        IF (NINCX.EQ.1) MRSURF=MRSURF+1
        GOTO 204
      ENDIF
C
C  PARTICLE OUTSIDE STANDARD MESH
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
C****ELLIPTICAL MESH  (NOT NECESSAIRLY CONCENTRIC OR CONFOCAL)
C
C  COEFFICIENTS OF QUADRATIC
      IF (NJUMP.EQ.0) THEN
        IF (NLSRFX) THEN
          X0SURF=X0-EP1(MRSURF)
          DSRF=ELLQ(MRSURF)
          ZR=RQ(MRSURF)
          Y0Q=(ZR-X0SURF*X0SURF)*DSRF
          MSAVE=MRSURF
        ELSE
          MSAVE=0
          Y0Q=Y0*Y0
        ENDIF
        YVY=Y0*VELY
        VELXQ=VELX*VELX
        VELYQ=VELY*VELY
        XX0=X0
        VVELX=VELX
C   TEST FOR DIRECTION, TENTATIVELY ASSUME...
        IF (NRCELL.EQ.0) THEN
          WRITE (6,*) 'NRCELL=0 IN TIMER: NOT READY '
          CALL EXIT
        ENDIF
        NINCX=-1
        MRSURF=NRCELL
        ZB=(XX0-EP1(MRSURF))*VVELX+YVY/ELLQ(MRSURF)
        IF(ZB.LE.0.) GOTO 1200
        NINCX=1
        MRSURF=MRSURF+1
        GOTO 1100
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
C     IF (NLTRC) WRITE (6,*) 'IPOLGN FROM TIMER ',IPOLGN,NRCELL
      IF (NRI.EQ.0) THEN
        RETURN
      ELSE
        IF (PT.GT.0.AND.PT.LT.PS) THEN
          IS=MRSURF-MAX(0,NINCX)
          MS=MRSURF
          PS=PT
        ENDIF
        NRI=IR+1
        GOTO 2000
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
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(5X,A)')
     .    'TIMER: Input  (NRCELL,NPCELL,NJUMP,NLSRFX,NLSRFY,'//
     .    'NLSRFT MRSURF,MPSURF,IPOLG,NINCX) '
        WRITE(6,'(5X,3I4,2X,3L2,4I5)')
     .    NRCELL,NPCELL,NJUMP,NLSRFX,NLSRFY,NLSRFT,MRSURF,MPSURF,
     .    IPOLG,NINCX
      ENDIF
c slmod end
      IF (NRCELL.NE.0) THEN
c slmod begin - tr
        AP=0
c slmod end
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
c slmod begin - new
          IF (NINCX.GE. 1) LNGB1=.FALSE.
          IF (NINCX.LE.-1) LNGB3=.FALSE.

c...DEV
          IF (NINCX.EQ.-17) THEN
            LNGB1=.FALSE.
            LNGB3=.TRUE.
          ENDIF
          IF (NINCX.EQ.+17) THEN
            LNGB1=.TRUE.
            LNGB3=.FALSE.
          ENDIF
c
c          IF (NINCX.EQ.1) LNGB1=.FALSE.
c          IF (NINCX.EQ.-1) LNGB3=.FALSE.
c slmod end
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
c slmod begin - tr
c...BUG: Not sure what is going on here, but it seems as if
c        the particle intersects the target surface after
c        being launched.  The really strange thing is that 
c        it is still closer to the next poloidal surface than the 
c        target surface, which doesn't make much sense.  However,
c        this only affects a few particles, so I am not
c        spending much time on it here (although, I think something
c        similar may have been happening to an ASDEX grid that Karl
c        was using when David was there, where EIRENE seemed
c        to be launching a particle from behind the target).  Anyway,
c        I am not sure where the problem is, in DIVIMP or EIRENE, but
c        I am putting in this rather specific fix for now:
c TURNED OFF FOR MAST GRID...
        ELSEIF (.FALSE..AND.
     .          NLSRFT.AND.NPCELL.EQ.1.AND.NRCELL.EQ.11) THEN
          LNGB4=.FALSE.
          WRITE(6,*) 'WARNING: CRUDE FIX TO TARGET CELL PROBLEM',NPANU
c slmod end
        ENDIF
C
        NCOUP=0
C  CALCULATE INTERSECTIONS OF FLIGHT WITH CELL BOUNDARIES
6001    CONTINUE
c slmod begin - tr - new
        ILOOP=ILOOP+1

        IF (IP.EQ.0) THEN
          IF (output) WRITE(0,*) 'ERROR: IP = 0'
          WRITE(6,*) 'ERROR: IP = 0'
          WRITE(6,*) '       NPANU = ',npanu
          PT = 1.D30
          RETURN
        ENDIF

        IF (ILOOP.GT.N2ND) THEN
          IF (output) 
     .      WRITE(0,*) 'ERROR: LIKELY AN INFINITE LOOP IN TIMER'
          WRITE(6,*) 'ERROR: LIKELY AN INFINITE LOOP IN TIMER'
          WRITE(6,*) '       NPANU = ',npanu
          PT=1.D30
          RETURN
        ENDIF

        IF (RVRTAG(IR,IP).EQ.1) THEN
          LNGB1 = .FALSE.
          LNGB3 = .FALSE.
        ENDIF
c        IF (PVRTAG(IR,IP).EQ.1) THEN
c          LNGB2 = .FALSE.
c          LNGB4 = .FALSE.
c        ENDIF

c       IPP/08 Krieger - ensure index of pvrtag is not zero
        IF (NGHPOL(2,IR,IP).GT.0.AND.
     .      PVRTAG(IR,max(1,NGHPOL(2,IR,IP))).EQ.1) LNGB2 = .FALSE.

c        IF (PVRTAG(IR,NGHPOL(4,IR,IP)).EQ.1) LNGB4 = .FALSE.
c        IF (PVRTAG(IR,IP+1).EQ.1) LNGB2 = .FALSE.

c        IF (npanu.EQ.6) WRITE(6,*) 'DATA:',ir,ip,PVRTAG(IR,IP),lngb4

        IF (PVRTAG(IR,IP).EQ.1) LNGB4 = .FALSE.

c        IF (npanu.EQ.6) WRITE(6,*) 'DATA:',ir,ip,PVRTAG(IR,IP),lngb4


        IF (GRIDOPT.EQ.1) THEN
c...note: not sure what this IF block was supposed to be for
        ENDIF
c slmod end
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
c slmod begin - grid - tr
        IF (GRIDOPT.EQ.1) THEN
          IF (LNGB1) THEN
            VPLX1=XVERT(IR,IP+1,1)-XVERT(IR,IP,1)
            VPLY1=YVERT(IR,IP+1,1)-YVERT(IR,IP,1)
            T1=((XVERT(IR,IP,1)-X0)*VELY-
     .          (YVERT(IR,IP,1)-Y0)*VELX)/
     .         (VELX*VPLY1-VELY*VPLX1+EPS60)
          ENDIF
          IF (LNGB3) THEN
            VPLX3=XVERT(IR,IP+1,2)-XVERT(IR,IP,2)
            VPLY3=YVERT(IR,IP+1,2)-YVERT(IR,IP,2)
            T3=((XVERT(IR,IP,2)-X0)*VELY-
     .          (YVERT(IR,IP,2)-Y0)*VELX)/
     .         (VELX*VPLY3-VELY*VPLX3+EPS60)
          ENDIF
          IF (LNGB2)
     .      T2=((XVERT(IR,IP+1,1)-X0)*VELY-
     .          (YVERT(IR,IP+1,1)-Y0)*VELX)/
     .         (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
          IF (LNGB4)
     .      T4=((XVERT(IR,IP,1)-X0)*VELY-
     .          (YVERT(IR,IP,1)-Y0)*VELX)/
     .         (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)


c       WRITE(6,'(2X,4(3I4,2F10.6))')
c     .   ir,ip  ,1,XVERT(IR,IP  ,1),YVERT(IR,IP  ,1),
c     .   ir,ip+1,1,XVERT(IR,IP+1,1),YVERT(IR,IP+1,1),
c     .   ir,ip  ,2,XVERT(IR,IP  ,2),YVERT(IR,IP  ,2),
c     .   ir,ip+1,2,XVERT(IR,IP+1,2),YVERT(IR,IP+1,2)


c        WRITE(6,*) 'VELX,Y=',velx,vely,velz

c
c Debug:
c
c...note: leave this for now
          IF (LNGB1) THEN
            DUM=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .         (VELX*VPLY(IR,IP)-VELY*VPLX(IR,IP)+EPS60)
            CALL CHECKNUM('ZAKA1',IR,IP,T1,DUM)
          ENDIF

          IF (LNGB2) THEN
            DUM=((XPOL(IR,IP+1)-X0)*VELY-(YPOL(IR,IP+1)-Y0)*VELX)/
     .          (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
            CALL CHECKNUM('ZAKA2',IR,IP,T2,DUM)
          ENDIF

          IF (LNGB3) THEN
            DUM=((XPOL(IR+1,IP)-X0)*VELY-(YPOL(IR+1,IP)-Y0)*VELX)/
     .          (VELX*VPLY(IR+1,IP)-VELY*VPLX(IR+1,IP)+EPS60)
            CALL CHECKNUM('ZAKA3',IR,IP,T3,DUM)
          ENDIF

          IF (LNGB4) THEN
            DUM=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .          (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)
            CALL CHECKNUM('ZAKA4',IR,IP,T4,DUM)
          ENDIF



        ELSE
          IF (LNGB1)
     .    T1=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .       (VELX*VPLY(IR,IP)-VELY*VPLX(IR,IP)+EPS60)
          IF (LNGB2)
     .    T2=((XPOL(IR,IP+1)-X0)*VELY-(YPOL(IR,IP+1)-Y0)*VELX)/
     .        (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
          IF (LNGB3)
     .    T3=((XPOL(IR+1,IP)-X0)*VELY-(YPOL(IR+1,IP)-Y0)*VELX)/
     .       (VELX*VPLY(IR+1,IP)-VELY*VPLX(IR+1,IP)+EPS60)
          IF (LNGB4)
     .    T4=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .        (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)
        ENDIF
c
c        IF (LNGB1)
c     .  T1=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
c     .     (VELX*VPLY(IR,IP)-VELY*VPLX(IR,IP)+EPS60)
c        IF (LNGB2)
c     .  T2=((XPOL(IR,IP+1)-X0)*VELY-(YPOL(IR,IP+1)-Y0)*VELX)/
c     .      (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
c        IF (LNGB3)
c     .  T3=((XPOL(IR+1,IP)-X0)*VELY-(YPOL(IR+1,IP)-Y0)*VELX)/
c     .     (VELX*VPLY(IR+1,IP)-VELY*VPLX(IR+1,IP)+EPS60)
c        IF (LNGB4)
c     .  T4=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
c     .      (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)
c slmod end
C       IF (NLTRC) WRITE (6,*) ' T1,T2,T3,T4 ',T1,T2,T3,T4
        LCT1 = T1.GE.0.D0 .AND. T1.LE.1.D0
        LCT2 = T2.GE.0.D0 .AND. T2.LE.1.D0
        LCT3 = T3.GE.0.D0 .AND. T3.LE.1.D0
        LCT4 = T4.GE.0.D0 .AND. T4.LE.1.D0
C  CALCULATE TIME OF FLIGHT FROM STARTING POINT X=(X0,Y0,Z0)
C  TO THE BOUNDARY OF THE ACTUELL CELL
c slmod begin - tr
        IF (GRIDOPT.EQ.1) THEN
          VPLX1=XVERT(IR,IP+1,1)-XVERT(IR,IP,1)
          VPLY1=YVERT(IR,IP+1,1)-YVERT(IR,IP,1)
          VPLX3=XVERT(IR,IP+1,2)-XVERT(IR,IP,2)
          VPLY3=YVERT(IR,IP+1,2)-YVERT(IR,IP,2)
          VVTX2=VVTX(IR,IP+1)
          VVTY2=VVTY(IR,IP+1)
          VVTX4=VVTX(IR,IP)
          VVTY4=VVTY(IR,IP)
          IF (ABS(VELX).GT.ABS(VELY)) THEN
            IF (LCT1) PT1=(XVERT(IR,IP  ,1)-X0+VPLX1*T1)/VELX
            IF (LCT2) PT2=(XVERT(IR,IP+1,1)-X0+VVTX2*T2)/VELX
            IF (LCT3) PT3=(XVERT(IR,IP  ,2)-X0+VPLX3*T3)/VELX
            IF (LCT4) PT4=(XVERT(IR,IP  ,1)-X0+VVTX4*T4)/VELX
          ELSE
            IF (LCT1) PT1=(YVERT(IR,IP  ,1)-Y0+VPLY1*T1)/VELY
            IF (LCT2) PT2=(YVERT(IR,IP+1,1)-Y0+VVTY2*T2)/VELY
            IF (LCT3) PT3=(YVERT(IR,IP  ,2)-Y0+VPLY3*T3)/VELY
            IF (LCT4) PT4=(YVERT(IR,IP  ,1)-Y0+VVTY4*T4)/VELY
          ENDIF

c          IF (ir.EQ.23.AND.ip.EQ.17) THEN
c            WRITE(6,*) '???',LCT1,LCT2,LCT3,LCT4
c            WRITE(6,*) '???',ABS(VELX),ABS(VELY)
c            WRITE(6,*) '???',XVERT(IR,IP,1),X0,VPLX1,VELX
c          ENDIF

c
c Debug:
c
          IF (ABS(VELX).GT.ABS(VELY)) THEN
            IF (LCT1) THEN
              DUM=(XPOL(IR,IP)-X0+VPLX(IR,IP)*T1)/VELX
              CALL CHECKNUM('ZAQA1',IR,IP,PT1,DUM)
            ENDIF
            IF (LCT2) THEN
              DUM=(XPOL(IR,IP+1)-X0+VVTX(IR,IP+1)*T2)/VELX
              CALL CHECKNUM('ZAQA2',IR,IP,PT2,DUM)
            ENDIF
            IF (LCT3) THEN
              DUM=(XPOL(IR+1,IP)-X0+VPLX(IR+1,IP)*T3)/VELX
              CALL CHECKNUM('ZAQA3',IR,IP,PT3,DUM)
            ENDIF
            IF (LCT4) THEN
              DUM=(XPOL(IR,IP)-X0+VVTX(IR,IP)*T4)/VELX
              CALL CHECKNUM('ZAQA4',IR,IP,PT4,DUM)
            ENDIF
          ELSE
            IF (LCT1) THEN
              DUM=(YPOL(IR,IP)-Y0+VPLY(IR,IP)*T1)/VELY
              CALL CHECKNUM('ZAQA5',IR,IP,PT1,DUM)
            ENDIF
            IF (LCT2) THEN
              DUM=(YPOL(IR,IP+1)-Y0+VVTY(IR,IP+1)*T2)/VELY
              CALL CHECKNUM('ZAQA6',IR,IP,PT2,DUM)
            ENDIF
            IF (LCT3) THEN
              DUM=(YPOL(IR+1,IP)-Y0+VPLY(IR+1,IP)*T3)/VELY
              CALL CHECKNUM('ZAQA7',IR,IP,PT3,DUM)
            ENDIF
            IF (LCT4) THEN
              DUM=(YPOL(IR,IP)-Y0+VVTY(IR,IP)*T4)/VELY
              CALL CHECKNUM('ZAQA8',IR,IP,PT4,DUM)
            ENDIF
          ENDIF

        ELSE
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
        ENDIF
c
c        IF (ABS(VELX).GT.ABS(VELY)) THEN
c          IF (LCT1) PT1=(XPOL(IR,IP)-X0+VPLX(IR,IP)*T1)/VELX
c          IF (LCT2) PT2=(XPOL(IR,IP+1)-X0+VVTX(IR,IP+1)*T2)/VELX
c          IF (LCT3) PT3=(XPOL(IR+1,IP)-X0+VPLX(IR+1,IP)*T3)/VELX
c          IF (LCT4) PT4=(XPOL(IR,IP)-X0+VVTX(IR,IP)*T4)/VELX
c        ELSE
c          IF (LCT1) PT1=(YPOL(IR,IP)-Y0+VPLY(IR,IP)*T1)/VELY
c          IF (LCT2) PT2=(YPOL(IR,IP+1)-Y0+VVTY(IR,IP+1)*T2)/VELY
c          IF (LCT3) PT3=(YPOL(IR+1,IP)-Y0+VPLY(IR+1,IP)*T3)/VELY
c          IF (LCT4) PT4=(YPOL(IR,IP)-Y0+VVTY(IR,IP)*T4)/VELY
c        ENDIF
c slmod end
        LCT1 = LCT1 .AND. PT1.GE.0.D0
        LCT2 = LCT2 .AND. PT2.GE.0.D0
        LCT3 = LCT3 .AND. PT3.GE.0.D0
        LCT4 = LCT4 .AND. PT4.GE.0.D0
c slmod begin - debug - tr
        IF (.NOT.(LCT1.OR.LCT2.OR.LCT3.OR.LCT4))THEN
          IF (T1.GE.-1.0D-05.AND.T1.LE.0.0D0) THEN
            IF (output) WRITE(0,*) 'ALMOST 1: NPANU =',npanu
            WRITE(6,*) 'ALMOST 1: NPANU =',npanu
          ENDIF
          IF (T2.GE.-1.0D-05.AND.T2.LE.0.0D0) THEN
            IF (output) WRITE(0,*) 'ALMOST 2: NPANU =',npanu
            WRITE(6,*) 'ALMOST 2: NPANU =',npanu
          ENDIF
          IF (T3.GE.-1.0D-05.AND.T3.LE.0.0D0) THEN
            IF (output) WRITE(0,*) 'ALMOST 3: NPANU =',npanu
            WRITE(6,*) 'ALMOST 3: NPANU =',npanu
          ENDIF
          IF (T4.GE.-1.0D-05.AND.T4.LE.0.0D0) THEN
            IF (output) WRITE(0,*) 'ALMOST 4: NPANU =',npanu
            WRITE(6,*) 'ALMOST 4: NPANU =',npanu
          ENDIF
        ENDIF

        IF (PRINTOPT.EQ.1) THEN

          WRITE(6,'(5X,A,2I4,I2,1X,4L2,1X,4L2,1X,1P,4E10.2,
     .              2X,4E10.2)')
     .      'TIMER: ',
     .      ir,ip,rvrtag(ir,ip),lngb1,lngb2,lngb3,lngb4,
     .      lct1,lct2,lct3,lct4,t1,t2,t3,t4,pt1,pt2,pt3,pt4
        ENDIF
c slmod end
C       IF (NLTRC) WRITE (6,*) ' LCT1,LCT2,LCT3,LCT4 ',
C    .                           LCT1,LCT2,LCT3,LCT4
C       IF (NLTRC) WRITE (6,*) ' PT1,PT2,PT3,PT4 ',PT1,PT2,PT3,PT4
        T=MAX(PT1,PT2,PT3,PT4)
C       IF (NLTRC) WRITE (6,*) ' T = ',T
C  IF INTERSECTION WITH POLOIDAL BOUNDARY CONTINUE WITH NEIGHBORING CELL
c slmod begin - tr
          IF (LCT2.AND.LCT4.AND.LCT1) THEN
            IF (output) WRITE(0,*) 'DISCARDING INTERSECTIONS A'
            WRITE(6,*) 'DISCARDING INTERSECTIONS A'
            LCT2 = .FALSE.
            LCT4 = .FALSE.
            T=MAX(PT1,PT3)
          ENDIF
          IF (LCT2.AND.LCT4.AND.LCT3) THEN
            IF (output) WRITE(0,*) 'DISCARDING INTERSECTIONS B'
            WRITE(6,*) 'DISCARDING INTERSECTIONS A'
            LCT2 = .FALSE.
            LCT4 = .FALSE.
            T=MAX(PT1,PT3)
          ENDIF
c slmod end
        IF (LCT2.OR.LCT4) THEN
c slmod begin - tr

          IF (LCT2.AND.LCT4) THEN
c          IF (LCT2.AND.LCT4.AND..NOT.NLSRFT) THEN
            IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10) THEN
              WRITE(6,'(5X,A)') 'TIMER: TROUBLE - 2 POLOIDAL '//
     .                          'INTERSECTIONS FOUND'
            ENDIF

c              isl1 = -1
c              isl2 = -1
c              isl1 = LEARC1(x00,y00,0.0D0,isl2,nrcell,nrcell,
c     .                      .TRUE.,.FALSE.,
c     .                      npanu,'TEST 03     ')
c              WRITE(6,'(5X,A,I5)') 'TIMER: LEARC1 = ',isl2
c

c...kill particle
            IF (output) 
     .        WRITE(0,*) 'TIMER: 2 INTERSECTIONS, KILLING PARTICLE'
            WRITE(6,*) 'TIMER: 2 INTERSECTIONS, KILLING PARTICLE'
            WRITE(6,*) '       NPANU = ',npanu
            PT=1.D30
            RETURN

          ENDIF
c slmod end
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
          IF (IP.NE.0) GOTO 6001
C  NO NEIGHBORING CELL: PARTICLE HAS HIT A POLOIDAL BOUNDARY OF THE MESH
c slmod begin - debug - tr
          IF (PRINTOPT.EQ.1)
     .      WRITE(6,'(5X,A)') 'TIMER: PARTICLE HAS HIT A POLOIDAL '//
     .                        'BOUNDARY OF THE MESH'
c slmod end
          MRSURF=0
          PT=1.D30
          NINCX=0
        ELSEIF (LCT1.OR.LCT3) THEN
C  INTERSECTION WITH RADIAL CELL BOUNDARY FOUND
c slmod begin - debug - tr
          IF (PRINTOPT.EQ.1)
     .      WRITE(6,'(5X,A)') 'TIMER: INTERSECTION WITH RADIAL '//
     .                        'CELL BOUNDARY FOUND'
c slmod end
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
c slmod begin - new
c...DEV
          IF (OPTCONMAP.EQ.1) THEN
            IF (LCT1) THEN
              MRSURF=IR
              NINCX=NGHPOL(1,IR,IP)-IR
              ICELLR=NGHPOL(1,IR,IP)
              IF (IR.GT.1.AND.ICELLR.LE.0) THEN
c                WRITE(0,*) 'TROUBLE DATA A:',IR,IP,NGHPOL(1,IR,IP),NINCX
c                WRITE(6,*) 'TROUBLE DATA A:',IR,IP,NGHPOL(1,IR,IP),NINCX
c                STOP 'TROUBLE! A'
              ENDIF
c            WRITE(0,*) 'LOWER: ',nincx,-1,icellr,ir-1
            ELSE
c...          This could be a problem, and it seems to me that it should have 
c             been for the PFZ extension on the C-Mod grid. I must have put in a
c             fix somewhere:
              MRSURF=NGHPOL(3,IR,IP)
              ICELLR=NGHPOL(3,IR,IP)
              NINCX=NGHPOL(3,IR,IP)-IR
              IF (ICELLR.LE.0) THEN
                IF (IR.EQ.4) THEN  ! HACK, CHEAT, TEMP, FOR CMOD GRID 
                  MRSURF=5
                  ICELLR=IP
                  NINCX=+1
                  IF (funnyoutput) THEN
                    WRITE(0,*)
                    WRITE(0,*) '*************************************'
                    WRITE(0,*) '* C-MOD GRID CUSTOM OVERRIDE IN USE *'
                    WRITE(0,*) '*************************************'
                    funnyoutput = .FALSE.
                  ENDIF
c                  WRITE(6,*) 'FUNNY ONE'
                ELSE
                  WRITE(0,*) 'TROUBLE DATA B:',IR,IP,NGHPOL(3,IR,IP),NINCX
                  WRITE(6,*) 'TROUBLE DATA B:',IR,IP,NGHPOL(3,IR,IP),NINCX
                  STOP 'TROUBLE! B'
                ENDIF
              ENDIF
c              WRITE(0,*) 'HIGHER:',nincx,+1,icellr,ir+1
            ENDIF
          ELSE
            IF (LCT1) THEN
              MRSURF=IR
              NINCX=-1
              ICELLR=IR-1
            ELSE
              MRSURF=IR+1
              ICELLR=IR+1
              NINCX=1
            ENDIF
          ENDIF
c
c          IF (LCT1) THEN
c            MRSURF=IR
c            NINCX=-1
c            ICELLR=IR-1
c          ELSE
c            MRSURF=IR+1
c            ICELLR=IR+1
c            NINCX=1
c          ENDIF
c slmod end
        ELSE
C  NO INTERSECTION FOUND
c slmod begin - tr - new
          IP2=NGHPOL(2,IR,IP)
          IF (((AP.EQ.0.OR.AP.EQ.-1).AND.PVRTAG(IR,IP ).EQ.1).OR.
     .        ((AP.EQ.0.OR.AP.EQ.+1).AND.PVRTAG(IR,IP2).EQ.1)) THEN

            LNGB1=.TRUE.
            LNGB2=.TRUE.
            LNGB3=.TRUE.
            LNGB4=.TRUE.

c            IF     (AP.EQ.0.AND.IP.EQ.NPOINT(1,1)) THEN
c              AP=+1
c              LNGB4=.FALSE.
c              IP=NGHPOL(2,IR,IP)
c            ELSEIF (AP.EQ.0.AND.IP.EQ.NPOINT(2,NPPLG)) THEN
c              AP=-1
c              LNGB4=.FALSE.
c              IP=NGHPOL(2,IR,IP)
            IF     (NGHPOL(2,IR,IP).GT.0.AND.
     .              PVRTAG(IR,NGHPOL(2,IR,IP)).EQ.1.AND.AP.EQ.0) THEN
              AP=+1
            ELSEIF (NGHPOL(4,IR,IP).GT.0.AND.
     .              PVRTAG(IR,NGHPOL(4,IR,IP)).EQ.1.AND.AP.EQ.0) THEN
              AP=-1
            ENDIF

            IF (AP.EQ.+1) THEN
c              WRITE(6,*) 'MOVING A:',IP,NGHPOL(2,IR,IP)
              IP=NGHPOL(2,IR,IP)
              LNGB4=.FALSE.
            ELSE
c              WRITE(6,*) 'MOVING B:',IP,NGHPOL(4,IR,IP)
              IP=NGHPOL(4,IR,IP)
              LNGB2=.FALSE.
            ENDIF

c            IF (NLSRFY) THEN
c              IF     (IP    .EQ.NPOINT(1,1).AND.
c     .                MPSURF.EQ.NPOINT(1,1)) THEN
c                AP=+1
c                LNGB4=.FALSE.
c              ELSEIF (IP    .EQ.NPOINT(2,NPPLG).AND.
c     .                MPSURF.EQ.NPOINT(2,NPPLG)) THEN
c                AP=-1
c                LNGB2=.FALSE.
c              ENDIF
c            ELSEIF (IP.EQ.NPOINT(1,1)) THEN
c              WRITE(0,*) 'MOVING E:',ip,NGHPOL(2,IR,IP)
c              AP=+1
c            ELSEIF (IP.EQ.NPOINT(2,NPPLG)) THEN
c              WRITE(0,*) 'MOVING F:',ip,NGHPOL(4,IR,IP)
c              AP=-1
c            ENDIF

c MOVE POLOIDAL CELL INDEX
c            IF (AP.EQ.+1) THEN
c              WRITE(6,*) 'MOVING A:',IP,NGHPOL(2,IR,IP)
c              IP=NGHPOL(2,IR,IP)
c            ELSE
c              WRITE(6,*) 'MOVING B:',IP,NGHPOL(4,IR,IP)
c              IP=NGHPOL(4,IR,IP)
c            ENDIF


            NJUMP=1
            GOTO 6001

c            IF (NLSRFY.AND.IP.EQ.NPOINT(2,1)-1) THEN
c              IF (MPSURF.EQ.NPOINT(2,1)) THEN
c                IF (ir.EQ.23) WRITE(0,*) 'MOVING A:',ip,NGHPOL(2,IR,IP)
c                IP=NGHPOL(2,IR,IP)
c                LNGB4=.FALSE.
c                AP=+1
c              ELSE
c                IF (ir.EQ.23) WRITE(0,*) 'MOVING B:',ip,NGHPOL(4,IR,IP)
c     .          ,mpsurf
c                IP=NGHPOL(4,IR,IP)
c                LNGB2=.FALSE.
c                AP=-1
c              ENDIF
c            ELSEIF (NLSRFY.AND.IP.EQ.NPOINT(1,NPPLG)) THEN
c              IF (MPSURF.EQ.NPOINT(1,NPPLG)) THEN
c                IF (ir.EQ.23) WRITE(0,*) 'MOVING C',ip,NGHPOL(4,IR,IP)
c                IP=NGHPOL(4,IR,IP)
c                LNGB2=.FALSE.
c                AP=-1
c              ELSE
c                IF (ir.EQ.23) WRITE(0,*) 'MOVING D',ip,NGHPOL(2,IR,IP)
c                IP=NGHPOL(2,IR,IP)
c                LNGB4=.FALSE.
c                AP=+1
c              ENDIF
c            ELSEIF (PVRTAG(IR,IP).EQ.0) THEN
c              IF (ir.EQ.23) WRITE(0,*) 'MOVING E:',ip,NGHPOL(2,IR,IP)
c              IP=NGHPOL(2,IR,IP)
c              AP=+1
c            ELSEIF (PVRTAG(IR,IP+1).EQ.0) THEN
c              IF (ir.EQ.23) WRITE(0,*) 'MOVING F:',ip,NGHPOL(4,IR,IP)
c              IP=NGHPOL(4,IR,IP)
c              AP=-1
c           ELSE
c              IF (ir.EQ.23) WRITE(0,*) 'MOVING G:',ip,NGHPOL(4,IR,IP)
c              IP=IP+AP
c            ENDIF
c            NJUMP=1
c            GOTO 6001
          ENDIF

          IF (PRINTOPT.EQ.1)
     .      WRITE(6,'(5X,A)') 'TIMER: NO INTERSECTION FOUND '
c slmod end
          IF (NLSRFY.AND.NJUMP.EQ.0) THEN
C  PLAY SAVE: TRY ONCE AGAIN, IF PARTICLE ON POL. SURFACE
            WRITE (6,*) ' NO INTERSECTION IN TIMER. TRY ONCE AGAIN '
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
c slmod begin - debug - tr (modified)
          GOTO 7999
c        RETURN
c slmod end
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
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10)
     .  WRITE(6,'(5X,A)') 'TIMER: CHECKING FOR RE-ENTRY'
c slmod end
      DO 6100 I=1,NR1ST
        ISTS=INMP1I(I,0,0)
        IF (PRINTOPT.EQ.1)
     .     WRITE(6,'(5X,A,2I6,E10.2)') 'TIMER: ISTS=',i,ists,timint(i)
C  SURFACE I IS NOT A NON DEFAULT RADIAL STANDARD SURFACE
        IF (ISTS.EQ.0) GOTO 6100
        IF (TIMINT(I).EQ.0) THEN
          TIMINT(I)=1
          NTIM(I)=0
c slmod begin - grid - tr
          IF (GRIDOPT.EQ.1) THEN
c...note: Arbitrary - fix
            IF (I.GT.DIVSUR.AND.
     .          .NOT.(OPTCONMAP.EQ.1.AND.RADMAP(I).EQ.0)) THEN
              IR1=I-1
              VP1=2
            ELSE
              IR1=I
              VP1=1
            ENDIF
          ENDIF
c slmod end
C
C   SEARCH FOR ALL POSSIBLE INTERSECTIONS WITHIN THE CELL
C
          DO 6011 J=1,NRPLG
6011        LCUT(J)=.FALSE.
C
          DO 6012 J=1,NPPLG
            DO 6012 K=NPOINT(1,J),NPOINT(2,J)-1
c slmod begin - grid - tr - new
c...note: hope this is okay!
              IF (INMP1I(I,K,0).EQ.0) GOTO 6012

              IF (GRIDOPT.EQ.1) THEN
c...debug:


                V1=(YVERT(IR1,K+1,VP1)-Y0)*VELX-
     .             (XVERT(IR1,K+1,VP1)-X0)*VELY
                V2=(YVERT(IR1,K  ,VP1)-Y0)*VELX-
     .             (XVERT(IR1,K  ,VP1)-X0)*VELY

                DUM=(YPOL(I,K+1)-Y0)*VELX-(XPOL(I,K+1)-X0)*VELY
                CALL CHECKNUM('ZEWA1',IR,IP,V1,DUM)

                DUM=(YPOL(I,K)-Y0)*VELX-(XPOL(I,K)-X0)*VELY
                CALL CHECKNUM('ZEWA2',IR,IP,V2,DUM)
              ELSE
                V1=(YPOL(I,K+1)-Y0)*VELX-(XPOL(I,K+1)-X0)*VELY
                V2=(YPOL(I,K)-Y0)*VELX-(XPOL(I,K)-X0)*VELY
              ENDIF
c
c              V1=(YPOL(I,K+1)-Y0)*VELX-(XPOL(I,K+1)-X0)*VELY
c              V2=(YPOL(I,K)-Y0)*VELX-(XPOL(I,K)-X0)*VELY
c slmod end
              LCUT(K)=V1*V2.LE.0.

                IF (printopt.EQ.1)
     .            WRITE(6,'(5X,A,5I4,4F10.4,3X,2F10.4,1P,4E10.2,0P,L4)')
     .              'TIMER: Searching (I,K,INMP1I) ',
     .              i,k,inmp1i(i,k,0),ir1,vp1,
     .              XVERT(IR1,K+1,VP1),XVERT(IR1,K  ,VP1),
     .              YVERT(IR1,K+1,VP1),YVERT(IR1,K  ,VP1),
     .              X0,Y0,VELX,VELY,V1,V2,lcut(k)



6012      CONTINUE
C
          IF (I.EQ.MPOL) LCUT(JPOL)=.FALSE.
          KAN=ILLZ(NRPLG,LCUT,1)+1
          KEN=NRPLG-ILLZ(NRPLG,LCUT,-1)
c sltmp
          IF (printopt.EQ.1)
     .      WRITE(6,'(5X,A,2I4)')
     .        'TIMER: KAN,KEN: ',
     .        kan,ken
C
C   COMPUTE THE FLIGHT TIMES TO THE INTERSECTION POINTS
C
          DO 6013 K=KAN,KEN
            IF (LCUT(K)) THEN
c slmod begin - tr
              IF (GRIDOPT.EQ.1) THEN
                VPLX1=XVERT(IR1,K+1,VP1)-XVERT(IR1,K,VP1)
                VPLY1=YVERT(IR1,K+1,VP1)-YVERT(IR1,K,VP1)
                T1=((XVERT(IR1,K,VP1)-X0)*VPLY1-
     .              (YVERT(IR1,K,VP1)-Y0)*VPLX1)/
     .             (VELX*VPLY1-VELY*VPLX1+EPS60)

                DUM=((XPOL(I,K)-X0)*VPLY(I,K)-(YPOL(I,K)-Y0)*VPLX(I,K))
     .              /(VELX*VPLY(I,K)-VELY*VPLX(I,K)+EPS60)
                CALL CHECKNUM('ZEWA4',I,K,T1,DUM)
              ELSE
                T1=((XPOL(I,K)-X0)*VPLY(I,K)-(YPOL(I,K)-Y0)*VPLX(I,K))
     .             /(VELX*VPLY(I,K)-VELY*VPLX(I,K)+EPS60)
              ENDIF
c
c              T1=((XPOL(I,K)-X0)*VPLY(I,K)-(YPOL(I,K)-Y0)*VPLX(I,K))
c     .           /(VELX*VPLY(I,K)-VELY*VPLX(I,K)+EPS60)
c slmod end
              IF (T1.GT.0.) THEN
                NTIM(I)=NTIM(I)+1
                TIMPOL(I,NTIM(I))=T1
                IIMPOL(I,NTIM(I))=K
              ENDIF
            ENDIF
c slmod begin - new
            IF (printopt.GE.1.AND.printopt.LE.10) THEN
              WRITE(6,'(5X,A,2I4,E12.5)')
     .          'TIMER: TIME',I,ISTS,T1
            ENDIF
c slmod end
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
        NINCX=SIGN(1.D0,VELX*PLNX(MRSURF,IPOLGN)+
     .        VELY*PLNY(MRSURF,IPOLGN))
      ENDIF
C
      NJUMP=1
      IPOLGO=IPOLGN
C
C     IF (NLTRC) WRITE (6,*) 'PT,MRSURF,IPOLGN,NINCX ',
C    .                        PT,MRSURF,IPOLGN,NINCX
c slmod begin - debug - tr (modified)
 7999 IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(5X,A,2I3,2I3,A,E12.5,2I3)')
     .    'TIMER: Output (NR,NJUMP,MR,IPOLGN,PT,NINCX,IR)',
     .    NRCELL,NJUMP,MRSURF,IPOLGN,'*',PT,NINCX,IRCELL

        WRITE(6,'(5X,A)') 'TIMER: '
      ENDIF
c slmod end
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
C  LOOP OVER ALL TRIANGLES ALONG TRACK, UNTIL FIRST BOUNDARY SURFACE
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
        NINCX=SIGN(1,INMTI(IPOLGN,NRCELL))
      ELSE
        GOTO 9999
      ENDIF
C     IF (NLTRC) WRITE (6,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
C    .                         NRCELL,MRSURF,PT,NINCX,IPOLGN

      RETURN
C
10000 CONTINUE
C
C  GENERAL GEOMETRY OPTION: PROVIDE FLIGHT TIME IN CURRENT CELL
C
      CALL TIMUSR(NRCELL,X0,Y0,Z0,VELX,VELY,VELZ,NJUMP,
     .            NEWCEL,TIM,ICOS,IERR)
      IF (IERR.NE.0) GOTO 9999
      PT=TIM
      IF (NEWCEL.GT.0) THEN
        NINCX=NEWCEL-NRCELL
        INMP1I(MRSURF,1,1)=0
      ELSEIF (NEWCEL.LT.0) THEN
        NINCX=ICOS
        MRSURF=NRCELL
        INMP1I(MRSURF,1,1)=-NEWCEL
      ELSEIF (NEWCEL.EQ.0) THEN
        WRITE (6,*) 'NEWCEL=0, EXIT FROM MESH '
        CALL EXIT
      ENDIF
      NLSRFX=.FALSE.
      RETURN
C
9999  CONTINUE
      WRITE (6,*) 'ERROR IN TIMER, EXIT CALLED AT 9999'
      WRITE (6,*) 'NPANU ', NPANU
      CALL EXIT
      END
C
      SUBROUTINE TIMEA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C   1 ST INTERSECTION OF THE RAY X+T*VX,Y+T*VY,Z+T*VZ WITH ONE OF THE
C   ADDITIONAL SURFACES, DEFINED BY 2.ND ORDER EQUATIONS
C   IT IS ALSO CHECKED, WHETHER THIS INTERSECTION TAKES PLACE INSIDE THE
C   SPECIFIED BOUNDARIES OF THOSE SURFACES
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'
      INCLUDE 'CLGIN'
      INCLUDE 'CGRID'
      INCLUDE 'CTRIG'
      INCLUDE 'CCONA'
      INCLUDE 'CLOGAU'
      INCLUDE 'CESTIM'
      INCLUDE 'COMUSR'
C
      DIMENSION XB(3),XC(3),XD(3),XE(3),XF(3),XG(3)
      LOGICAL RLBNOT(NLIM),LGJ,NLTRC,BITGET
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
        IF (output) WRITE(0,*) 'MARK: LOOP TO 97'
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
        IF (output) WRITE(0,*) 'MARK: 2 POINTS'
        IF (RLB(J).GE.2.6.AND.RLB(J).LT.2.9) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE X CO-ORDINATE
          IF (RLB(J).LT.2.75) RLB(J)=1.
          IF (RLB(J).GE.2.75) RLB(J)=1.5
          XLIMS1(J,1)=P1(1,J)
          XLIMS2(J,1)=P2(1,J)
          DZ3=(P1(3,J)-P2(3,J))
          DY2=(P1(2,J)-P2(2,J))
          A0LM(J)=DZ3*P1(2,J)-DY2*P1(3,J)
          A1LM(J)=0.
          A2LM(J)=-DZ3
          A3LM(J)=DY2
          YLIMS1(J,1)=MIN(P1(2,J),P2(2,J))
          YLIMS2(J,1)=MAX(P1(2,J),P2(2,J))
          ZLIMS1(J,1)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(J,1)=MAX(P1(3,J),P2(3,J))
          IF (P1(2,J).EQ.P2(2,J)) THEN
            YLIMS1(J,1)=YLIMS1(J,1)-0.1
            YLIMS2(J,1)=YLIMS2(J,1)+0.1
          ENDIF
          IF (P1(3,J).EQ.P2(3,J)) THEN
            ZLIMS1(J,1)=ZLIMS1(J,1)-0.1
            ZLIMS2(J,1)=ZLIMS2(J,1)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND X-INDEPENDENCE OF SURFACE EQUATION
          P3(1,J)=1.D55
          DL=SQRT(DY2**2+DZ3**2)
          DX=XLIMS2(J,1)-XLIMS1(J,1)
          IF (DX.GT.1.D20) DX=XDF
          SAREA(J)=DL*DX
          GOTO 90
        ELSEIF (RLB(J).GE.2.3.AND.RLB(J).LT.2.6) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE Y CO-ORDINATE
          IF (RLB(J).LT.2.45) RLB(J)=1.
          IF (RLB(J).GE.2.45) RLB(J)=1.5
          YLIMS1(J,1)=P1(2,J)
          YLIMS2(J,1)=P2(2,J)
          DZ3=(P1(3,J)-P2(3,J))
          DX1=(P1(1,J)-P2(1,J))
          A0LM(J)=DZ3*P1(1,J)-DX1*P1(3,J)
          A1LM(J)=-DZ3
          A2LM(J)=0.
          A3LM(J)=DX1
          XLIMS1(J,1)=MIN(P1(1,J),P2(1,J))
          XLIMS2(J,1)=MAX(P1(1,J),P2(1,J))
          ZLIMS1(J,1)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(J,1)=MAX(P1(3,J),P2(3,J))
          IF (P1(1,J).EQ.P2(1,J)) THEN
            XLIMS1(J,1)=XLIMS1(J,1)-0.1
            XLIMS2(J,1)=XLIMS2(J,1)+0.1
          ENDIF
          IF (P1(3,J).EQ.P2(3,J)) THEN
            ZLIMS1(J,1)=ZLIMS1(J,1)-0.1
            ZLIMS2(J,1)=ZLIMS2(J,1)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND Y-INDEPENDENCE OF SURFACE EQUATION
          P3(2,J)=1.D55
          DL=SQRT(DX1**2+DZ3**2)
          DY=YLIMS2(J,1)-YLIMS1(J,1)
          IF (DY.GT.1.E20) DY=YDF
          SAREA(J)=DL*DY
          GOTO 90
        ELSEIF (RLB(J).GE.2.0.AND.RLB(J).LT.2.3) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE Z CO-ORDINATE
          IF (RLB(J).LT.2.15) RLB(J)=1.
          IF (RLB(J).GE.2.15) RLB(J)=1.5
          ZLIMS1(J,1)=P1(3,J)
          ZLIMS2(J,1)=P2(3,J)
          DY2=(P1(2,J)-P2(2,J))
          DX1=(P1(1,J)-P2(1,J))
          A0LM(J)=DY2*P1(1,J)-DX1*P1(2,J)
          A1LM(J)=-DY2
          A2LM(J)=DX1
          A3LM(J)=0.
          XLIMS1(J,1)=MIN(P1(1,J),P2(1,J))
          XLIMS2(J,1)=MAX(P1(1,J),P2(1,J))
          YLIMS1(J,1)=MIN(P1(2,J),P2(2,J))
          YLIMS2(J,1)=MAX(P1(2,J),P2(2,J))
          IF (P1(1,J).EQ.P2(1,J)) THEN
            XLIMS1(J,1)=XLIMS1(J,1)-0.1
            XLIMS2(J,1)=XLIMS2(J,1)+0.1
          ENDIF
          IF (P1(2,J).EQ.P2(2,J)) THEN
            YLIMS1(J,1)=YLIMS1(J,1)-0.1
            YLIMS2(J,1)=YLIMS2(J,1)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND Z-INDEPENDENCE OF SURFACE EQUATION
          P3(3,J)=1.D55
          DL=SQRT(DX1**2+DY2**2)
          IF (NLTRZ) THEN
            DZ=ZLIMS2(J,1)-ZLIMS1(J,1)
            IF (DZ.GT.1.D20) DZ=ZDF
          ELSEIF (NLTRA) THEN
            XS=(P1(1,J)+P2(1,J))*0.5+RMTOR
            IF (ILTOR(J).GT.0) THEN
              IF (.NOT.NLTOR) STOP 'MARK: NLTRA DEVELOPMENT NEEDED'
              DZ=ZLIMS2(J,1)-ZLIMS1(J,1)
C             IF (DZ.GT.1.D20) ????
            ELSEIF (ILTOR(J).EQ.0) THEN
              DZ=ZLIMS2(J,1)-ZLIMS1(J,1)
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
        IF (output) WRITE(0,*) 'MARK: PREPARING ARRAYS'
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
        IF (output) WRITE(0,*) 'MARK: TIMEA0 WARNINGS'
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
          CALL EXIT
        ENDIF
        GOTO 94
C
98      WRITE (6,*) 'WARNING FROM TIMEA0, SURFACE NO J IS TURNED OF'
        WRITE (6,*) 'BECAUSE ILL DEFINED '
        WRITE (6,*) 'J= ',J
        IGJUM0(J)=1
        IF (NLIMPB >= NLIMPS) THEN
          DO 99 IAB=0,NLIMI
            IGJUM1(IAB,J)=1
99        CONTINUE
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
        IF (output) WRITE(0,*) 'MARK: GOING IN'
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
      RETURN
C
      ENTRY TIMEA1(MSURF,NCELL,NTCELL,XX,YY,ZZ,TMT,
     .             VXX,VYY,VZZ,VV,
c slmod begin - tr - new
     .             MASURF,XR,YR,ZR,SG,TL,NLTRC,NACELL1,NTRSEG)
c
c     .             MASURF,XR,YR,ZR,SG,TL,NLTRC)
c slmod end
      LM1=NLIMII(NCELL)
C
      LM2=NLIMIE(NCELL)
C  PARTICLE ON STANDARD SURFACE?
      IF (MSURF.GT.NLIM) MSURF=0
C  SAVE INITIAL CO-ORDINATES
      XS=XX
      YS=YY
      ZS=ZZ
      TS=TMT
      VXS=VXX
      VYS=VYY
      VZS=VZZ
      VVS=VV
      NNTCLS=NTCELL
C  WORKING CO-ORDINATES
      X=XX
      Y=YY
      Z=ZZ
      T=TMT
      VX=VXX
      VY=VYY
      VZ=VZZ
      V=VV
      NNTCL=NTCELL
c slmod begin - tr (modified)
      IF ( OPTZMOTION.EQ.1.OR.
     .    (OPTZMOTION.EQ.2.AND.NACELL1.EQ.0).OR.
     .    (OPTZMOTION.EQ.3.AND.(NACELL1.EQ.0.OR.
     .     (XX.LT.62.5.AND.YY.GT.-61.0))).OR.
     .    (OPTZMOTION.EQ.4.AND.NACELL1.GT.0.AND.
     .     (XX.GT.62.5.OR.YY.LT.-61.0)).OR.
     .    (OPTZMOTION.EQ.5.AND.NACELL1.GT.0)) THEN
        VZS=0.0D0
        VZ=0.0D0
      ENDIF
c slmod end
C
      NLLLI=MSURF
C     IF (NLTRC) THEN
C       WRITE (6,*) 'TIMEA ',X,Y,Z,T,VX,VY,VZ,V
C       WRITE (6,*) 'MSURF,NTCELL ',MSURF,NTCELL
C     ENDIF
C
      TADD=0.
C
1000  TMIN=1.D30
      TL=1.D30
      MASURF=0
c slmod begin - tr
c
c
c
c
c
c
c      DO IREG=1,1
      DO IREG=1,IGJUM4(NCELL,0)

        IF     (IGJUM4(NCELL,IREG).EQ. 0) THEN
        ELSEIF (IGJUM4(NCELL,IREG).EQ.-1) THEN
          LM1=1			
          LM2=SBGKI-1		
        ELSEIF (IGJUM4(NCELL,IREG).EQ.-2) THEN
          LM1=EBGKI+1
          LM2=HADDI
        ELSEIF (IGJUM4(NCELL,IREG).EQ.-3) THEN
          LM1=HADDI+1
          LM2=HSTDI
        ELSEIF (IGJUM4(NCELL,IREG).EQ.-4) THEN
c          WRITE(0,*) '--',hstdi+1,nlimi
          LM1=HSTDI+1
          LM2=NLIMI
        ELSE
          LM1=IGJUM4(NCELL,IREG)
          LM2=LM1
        ENDIF

c slmod end
C
C  LOOP OVER SURFACE NUMBER, DO 100
C
      DO 100 J=LM1,LM2
        IF (IGJUM0(J).NE.0) GOTO 100
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM1(NLLLI,J) .NE. 0) GOTO 100
        ELSE
c            IF (BITGET(IGJUM1,0,NLIMPS,NLLLI,J,NBITS))
c     .        WRITE(0,*) 'MARK: YEP 1',NCELL,J
          IF (BITGET(IGJUM1,0,NLIMPS,NLLLI,J,NBITS)) GOTO 100
        END IF
        IF (NCELL.LE.NOPTIM) THEN
c slmod begin - igjum3 - tr (already added to EIRENE02)
          IF (NLIMPB >= NLIMPS) THEN
            IF (IGJUM3(NCELL,J).NE.0) GOTO 100
          ELSE
            IF (BITGET(IGJUM3,0,NOPTIM,NCELL,J,NBITS)) GOTO 100
          ENDIF
c
c          IF (IGJUM3(NCELL,J).NE.0) GOTO 100
c slmod end
        ENDIF
C
C  FOR NLTRA OPTION ONLY:
C  X,Z,VX AND VZ ARE GIVEN IN TOROIDAL CELL NNTCL,
C  TRANSFORM CO-ORDINATES FOR THIS TRACK FROM LOCAL SYSTEM NNTCL
C  TO THE LOCAL SYSTEM ILTOR(J), IN WHICH SURFACE J IS GIVEN
C  IF (ILTOR(J).LE.0) THIS SURFACE HAS TOROIDAL SYMMETRY
C
        IF (NLTRA) THEN
c slmod begin - tr
          IF (.NOT.NLTOR) THEN
            IF (ILTOR(J).EQ.0.OR.ILTOR(J).EQ.NTRSEG) THEN
              NNTCL=NNTCLS
              X=XS
              Z=ZS
              VX=VXS
              VZ=VZS
            ELSE            
              GOTO 100
            ENDIF
          ELSEIF (ILTOR(J).GT.0.AND.ILTOR(J).NE.NNTCL) THEN
c
c          IF (ILTOR(J).GT.0.AND.ILTOR(J).NE.NNTCL) THEN
c slmod end
            CALL FZRTOR (X,Z,NNTCL,XXR,THET,NTNEW,.FALSE.,0)
            CALL FZRTRI (X,Z,ILTOR(J),XXR,THET,NTNEW)
            ROT=2.*(NNTCL-ILTOR(J))*ALPHA
            VSAVE=VX
            VX=COS(ROT)*VSAVE-SIN(ROT)*VZ
            VZ=SIN(ROT)*VSAVE+COS(ROT)*VZ
            NNTCL=ILTOR(J)
C  X,Z,VX AND VZ ARE NOW GIVEN IN CELL ILTOR(J). SO ARE THE COEFFICIENTS
C  OF SURFACE NO. J. FIND INTERSECTION IN THIS LOCAL SYSTEM
          ELSEIF (ILTOR(J).EQ.0) THEN
C  TOROIDALLY SYMMETRIC SURFACE, SURFACE COEFFICIENTS ARE THE SAME
C  IN EACH TOROIDAL CELL, THUS ESPECIALLY IN CELL NNTCL
            NNTCL=NNTCLS
            X=XS
            Z=ZS
            VX=VXS
            VZ=VZS
          ENDIF
        ENDIF
C
C  FIND INTERSECTION TIME TMX WITH BOUNDARY NO. J
C
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
        IF (TMA.LE.EPS12.OR.TMI.GT.TMIN) GOTO 100
        IF (RLB(J).LT.0.) GOTO 21
        IF (RLB(J).GT.0.) GOTO 31
C
C  RLB(J) .EQ. 0.  STATEMENT 11---20
C
11      IF (TMI.LE.EPS12) GOTO 12
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        TMX=TMI
        GOTO 70
C
12      IF (TMA.GT.TMIN) GOTO 100
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
22      IF (TMA.GT.TMIN) GOTO 100
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
          TST=ALIMS(J,I)+XLIMS(J,I)*XN+YLIMS(J,I)*YN+ZLIMS(J,I)*ZN
          LGJ=TST.LE.0.
          IF (LGJ.AND.I.LT.ILIN(J)) GOTO 27
        ENDIF
        IF (LGJ.AND.ISCN(J).GT.0) THEN
          I=0
29        I=1+I
          TST=ALIMS0(J,I)+
     .        XN*(XLIMS1(J,I)+XN*XLIMS2(J,I)+YN*XLIMS3(J,I))+
     .        YN*(YLIMS1(J,I)+YN*YLIMS2(J,I)+ZN*ZLIMS3(J,I))+
     .        ZN*(ZLIMS1(J,I)+ZN*ZLIMS2(J,I)+XN*YLIMS3(J,I))
          LGJ=TST.LE.0.
          IF (LGJ.AND.I.LT.ISCN(J)) GOTO 29
        ENDIF
        IF (LGJ) THEN
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
        LGJ=XN.LE.XLIMS2(J,1).AND.XN.GE.XLIMS1(J,1).AND.
     .      YN.LE.YLIMS2(J,1).AND.YN.GE.YLIMS1(J,1).AND.
     .      ZN.LE.ZLIMS2(J,1).AND.ZN.GE.ZLIMS1(J,1)
        IF (RLBNOT(J)) LGJ=.NOT.LGJ
        TMX=TMI
        IF (LGJ) GOTO 70
C
32      IF (TMA.GT.TMIN) GOTO 100
        WR=A2+A1*(TMA+TMA)
36      XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
38      IF (RLB(J).LT.2.) THEN
          LGJ=XN.LE.XLIMS2(J,1).AND.XN.GE.XLIMS1(J,1).AND.
     .        YN.LE.YLIMS2(J,1).AND.YN.GE.YLIMS1(J,1).AND.
     .        ZN.LE.ZLIMS2(J,1).AND.ZN.GE.ZLIMS1(J,1)
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
        TMX=TMA
        IF (LGJ) GOTO 70
        GOTO 100
C
C   A1*TMX+A2=0
C
40      TMA=-A2/A1
C       IF (NLTRC) WRITE (6,*) 'TIMEA AT 40, J,TMA,A1,A2 ',J,TMA,A1,A2
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
        WR=-A2
        IF (RLB(J)) 26,16,36
C
C   A2*TMX+A3=0
C
50      IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
C       IF (NLTRC) WRITE (6,*) 'TIMEA AT 50, J,TMA,A2,A3 ',J,TMA,A2,A3
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
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
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
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
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
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
        IF (TMA.LE.EPS12.OR.TMA.GT.TMIN) GOTO 100
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
        NNJ=NNTCL
        VXJ=VX
        VZJ=VZ
        MASURF=J
        SG=SIGN(1.D0,WR)
C       IF (NLTRC) THEN
C         WRITE (6,*) 'TIMEA, TL,XR,YR,ZR,MASURF,SG '
C         WRITE (6,*)         TL,XR,YR,ZR,MASURF,SG
C       ENDIF
C
C  LOOP OVER SURFACE-INDEX FINISHED
C
100   CONTINUE
c slmod begin - tr
c
c
c
c
c
c
       ENDDO
c slmod end
C
C **********************************************************************
C
c slmod begin - tr
C  IF NO INTERSECTION FOUND, RETURN
      IF (MASURF.EQ.0) THEN
        IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10)
     .    WRITE(6,'(4X,A)') 'TIMEA1: NO INTERSECTION FOUND'
        RETURN
      ENDIF
C  INTERSECTION AT SURFACE NO. MASURF
C  IF NOT TRANSPARENT, RETURN
      IF (ILIIN(MASURF).GT.0) THEN
        IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10)
     .    WRITE(6,'(4X,A)') 'TIMEA1: NON-TRANSPARENT SURFACE'
        RETURN
      ENDIF
C  IF TRANSPARENT BUT WRONG SIDE, RETURN
      IF (ILSIDE(MASURF)*SG.LT.0) THEN
        IF (PRINTOPT.GE.1.AND.PRINTOPT.LE.10)
     .    WRITE(6,'(4X,A)') 'TIMEA1: TRANSPARENT BUT WRONG SIDE'
        RETURN
      ENDIF
c
cC  IF NO INTERSECTION FOUND, RETURN
c      IF (MASURF.EQ.0) RETURN
cC  INTERSECTION AT SURFACE NO. MASURF
cC  IF NOT TRANSPARENT, RETURN
c      IF (ILIIN(MASURF).GT.0) RETURN
cC  IF TRANSPARENT BUT WRONG SIDE, RETURN
c      IF (ILSIDE(MASURF)*SG.LT.0) RETURN
c slmod end
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
        NNTCL=NNJ
        NNTCLS=NNJ
        VX=VXJ
        VXS=VXJ
        VZ=VZJ
        VZS=VZJ
        GOTO 1000
      ELSEIF (ILIIN(MASURF).LT.0) THEN
C  TRANSPARENT, BUT SWITCH AND/OR SURFACE TALLIES
c slmod begin - debug - tr
        IF (printopt.GE.1.AND.printopt.LE.10) THEN
          WRITE(6,'(4X,A)') 'TIMEA1: TRANSPARENT BUT SWITCH'
          WRITE(6,'(4X,A,I5,E12.5)')
     .      'TIMEA1: Output (MASURF,TL) ',
     .      masurf,tl
          WRITE(6,'(4X,A)') 'TIMEA1: '
        ENDIF
c slmod end
        RETURN
      ENDIF
C
      END
C
C
      SUBROUTINE TIMEP (ZRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C   CALCULATE TIME SEGMENTS IN Y- OR POLOIDAL MESH
C
C   INPUT:
C
C       NRCELL = RADIAL CELL NUMBER, IN WHICH THIS TRACK OF LENGTH ZRAD
C                IS PERFORMED
C              = 0, IF TRACK OUTSIDE RADIAL MESH
C                MRSURF .NE. 0 THEN
C                REENTRY FOUND AT RADIAL SURFACE MRSURF.
C                OTHERWISE: REENTRY AT NON DEFAULT
C                POLOIDAL SURFACES IS SEARCHED FOR IN THIS CALL
C       NTCELL =
C       ITCELL =
C       NPCELL = POLOIDAL CELL INDEX OF POINT X00,Y00,Z00
C                AT WHICH THIS TRACK OF LENGTH ZRAD STARTS
C  WARNING: IF NLSRFY, NPCELL MAY BE WRONG
C       IPCELL =
C       X00,Y00,Z00: STARTING POINT FOR THIS TRACK
C
C       NCOUT,BLPD,NCOUNT
C
C   OUTPUT:
C
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
C              = 0  IF TRACK COMPLETELY OUSIDE STANDARD MESH
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'CCONA'
      INCLUDE 'CGEOM'
      INCLUDE 'CGRID'
      INCLUDE 'CPOLYG'
      INCLUDE 'CLOGAU'
      INCLUDE 'CLGIN'
      INCLUDE 'COMSOU'
c slmod begin - not tr
c...note: Not sure why this is needed.
      INCLUDE 'COMUSR'
c slmod end
      LOGICAL LCUTY(N2NDPLG),LCUTX(N1ST)
      DIMENSION NCOUNS(N2ND+N3RD)
      SAVE
C
C     IF (NLTRC) THEN
C       CALL LEER(1)
C       IF (NRCELL.GT.0) THEN
C         WRITE (6,*) 'TIMEP FROM INSIDE, NPANU ', NPANU
C         WRITE (6,*) 'ZRAD,NRCELL,NPCELL '
C         WRITE (6,*) ZRAD,NRCELL,NPCELL
C       ELSE
C         WRITE (6,*) 'TIMEP FROM OUTSIDE, NPANU ', NPANU
C         WRITE (6,*) 'MRSURF,MTSURF,ZRAD '
C         WRITE (6,*)  MRSURF,MTSURF,ZRAD
C       ENDIF
C     ENDIF
C
      ICOUT=1
      IADD=0
      ZRADS=ZRAD
      IPOLGS=IPOLGN
C
      IF (NLTOR) THEN
C       NCOUT=NCOUT
        ZRAD=BLPD(1)
        NTCELL=NCOUNT(1)
        IF (NCOUT.GT.1) IPOLGN=0
C       IF (NLTRC) WRITE (6,*) 'WG. NLTOR: ZRAD,NTCELL ',ZRAD,NTCELL
      ELSE
        NCOUT=1
C       ZRAD=ZRAD
C       NTCELL=1
      ENDIF
C
10000 CONTINUE
C
      IF (LEVGEO.EQ.3) THEN
c slmod begin - tr
        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(5X,A,5I4,3X,1P,E12.5)')
     .      'TIMEP: Input  (MPSURF,NR,NP,IR,IPOLGN ZRAD) ',
     .      MPSURF,NRCELL,NPCELL,IRCELL,IPOLGN,ZRAD
c slmod end
C
        IF (NRCELL.GT.0) GOTO 10
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT POLOIDAL SURFACES
C
        NCOUP=0
        ZRD=ZRAD
c slmod begin - debug - tr
        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(5X,A)') 'TIMEP: CHECKING FOR RE-ENTRY'
c slmod end
        DO 3 ISTS=1,NSTSI
          MPTEST=INUMP(ISTS,2)
          IF (MPTEST.NE.0) THEN
C  TEST POLOIDAL SURFACE NO. MPTEST FOR REENTRY
c slmod begin - debug - not tr
c                IF (printopt.GE.1.AND.printopt.LE.10)
c     .            WRITE(6,'(5X,A,I4)')
c     .              'TIMEP: Poloidal surface (MPTEST) ',
c     .              mptest
c slmod end
            DO 2 IR=1,NR1STM
              I1=IR+1
c slmod begin - tr
              IF (GRIDOPT.EQ.1) THEN
                V1=(YVERT(IR,MPTEST,1)-Y00)*VELX-
     .             (XVERT(IR,MPTEST,1)-X00)*VELY
                V2=(YVERT(IR,MPTEST,2)-Y00)*VELX-
     .             (XVERT(IR,MPTEST,2)-X00)*VELY
              ELSE
                V1=(YPOL(IR,MPTEST)-Y00)*VELX-(XPOL(IR,MPTEST)-X00)*VELY
                V2=(YPOL(I1,MPTEST)-Y00)*VELX-(XPOL(I1,MPTEST)-X00)*VELY
              ENDIF
c
c              V1=(YPOL(IR,MPTEST)-Y00)*VELX-(XPOL(IR,MPTEST)-X00)*VELY
c              V2=(YPOL(I1,MPTEST)-Y00)*VELX-(XPOL(I1,MPTEST)-X00)*VELY
c slmod end
              LCUTX(IR)=V1*V2.LE.0.
2           CONTINUE
            DO 4 IR=1,NR1STM
              IF (LCUTX(IR)) THEN
c slmod begin - tr
                IF (GRIDOPT.EQ.1) THEN
                  T1=((XVERT(IR,MPTEST,1)-X00)*VVTY(IR,MPTEST)-
     .                (YVERT(IR,MPTEST,1)-Y00)*VVTX(IR,MPTEST))
     .               /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)

                  DUM=((XPOL(IR,MPTEST)-X00)*VVTY(IR,MPTEST)-
     .                (YPOL(IR,MPTEST)-Y00)*VVTX(IR,MPTEST))
     .               /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)
                  CALL CHECKNUM('ZAAZ1',IR,MPTEST,T1,DUM)
                ELSE
                  T1=((XPOL(IR,MPTEST)-X00)*VVTY(IR,MPTEST)-
     .                (YPOL(IR,MPTEST)-Y00)*VVTX(IR,MPTEST))
     .               /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)
                ENDIF

c...GERMANY:    DON'T REGISTER A COLLISION WITH A *VERY* SHORT
c               COLLISION TIME, SINCE IT IS LIKELY FOR
c               THE PARTICLES PRESENT LOCATION:
                IF (T1.LT.EPS12) THEN
c                  WRITE(0,*) 'IGNORING POLOIDAL SURFACE COLLISION',NPANU
                  T1 = -1.0D0
                ENDIF
c
c                T1=((XPOL(IR,MPTEST)-X00)*VVTY(IR,MPTEST)-
c     .              (YPOL(IR,MPTEST)-Y00)*VVTX(IR,MPTEST))
c     .             /(VELX*VVTY(IR,MPTEST)-VELY*VVTX(IR,MPTEST)+EPS60)
c
c slmod end
C               IF (NLTRC) WRITE (6,*) 'IR,MPTEST,T1 ',IR,MPTEST,T1
                IF (T1.LT.0..OR.T1.GE.ZRD) GOTO 4
C               IF (NLTRC) WRITE (6,*) 'VALID INTERSECTION AT T1= ',T1
c slmod begin - tr
                IF (printopt.GE.1.AND.printopt.LE.10)
     .            WRITE(6,'(5X,A,2I4,1P,E12.4)')
     .              'TIMEP: *** VALID INTERSECTION *** (IR,MPTEST T1) ',
     .              IR,MPTEST,T1
c slmod end
                NCOUP=1
                LUPC(NCOUP)=MPTEST
                IRSAVE=IR
                MUPC(NCOUP)=SIGN(1.D0,VELX*PPLNX(IR,MPTEST)+
     .                                VELY*PPLNY(IR,MPTEST))
                JUPC(NCOUP)=1
                ZRD=T1
              ENDIF
4           CONTINUE
          ENDIF
3       CONTINUE
c slmod begin - tr (modified)
        IF (printopt.GE.1.AND.printopt.LE.10)
     .    WRITE(6,'(5X,A,3I4)')
     .      'TIMEP: NCOUP,LUPC,MPSURF=',NCOUP,LUPC(MAX(1,NCOUP)),MPSURF

c GERMANY: GETTING RID OF THE SECOND REQUIREMENT SO THAT PARTICLES
c          CAN BOUNCE OF THE BACK OF THE OUTER TARGET SEQUENTIALLY.
c          I MAY BE ABLE TO GET AWAY WITH THIS BECAUSE OF THE THE 
c          NEW CONSTRAINT ON T1 ABOVE (T1>1.0E-14 REQUIRED):
        IF (NCOUP.GT.0) THEN
c
c        IF (NCOUP.GT.0.AND.LUPC(MAX(1,NCOUP)).NE.MPSURF) THEN
c slmod end
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
C         IF (NLTRC) THEN
C           WRITE (6,*) 'REENTRY FOUND, MPSURF,ZRAD = ',MPSURF,ZRAD
C           WRITE (6,*) 'IRCELL ',IRCELL
C         ENDIF
c slmod begin - tr
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE (6,*) 'CHECK HERE FOR CORRECT POLOIDAL CELL NUMBER!'
c slmod end
          GOTO 31
        ELSE
C  NO REENTRY FOUND
          NCOUP=1
          JUPC(1)=1
          ALPD(1)=ZRAD
          IPCELL=IPOLGN
          MPSURF=0
C         IF (NLTRC) THEN
C           WRITE (6,*) 'NO REENTRY FOUND '
C         ENDIF
c slmod begin - debug - tr
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(5X,A,2F10.4)')
     .        'TIMEP: Updating position 5 (VELX,VELY) ',VELX,VELY
c slmod end
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
          WRITE (6,*) 'ERROR IN TIMEP: NCOUP=0'
c slmod begin - tr
c...note: Just trying to kill particle.
c          WRITE(0,*) 'TIMEP: KILL PARTICLE 01'
          WRITE(6,*) 'TIMEP: KILL PARTICLE  NPANU = ',NPANU
          ZRAD=1.0E+30
c slmod end
          RETURN
        ENDIF
C       IF (NLTRC) THEN
C         WRITE (6,*) ' TIMEP IN NEIGHBOR PART '
C         WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C         WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C         WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C         WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C       ENDIF
        NLSRFY=.FALSE.
        IF (SQRT((X0-X00)**2+(Y0-Y00)**2).GT.EPS10) THEN
          DO J=1,NCOUP
            ALPD(J)=ALPD(J)-ZT
          ENDDO
        ENDIF
C
C  ACCOUNT FOR ADDITIONAL SURFACES INSIDE MESH
C
        IF (ALPD(NCOUP).GT.ZRAD) THEN
          DO J=1,NCOUP
            IF (ALPD(J).GT.ZRAD) THEN
c slmod begin - tr (already there)
c... Bug fix sent by Detlev on March 16, 2001.
CPB
C  READJUST POLOIDAL SECTION IN CELL J
C  SAVE OTHER POLOIDAL SECTIONS ALPD FOR LATER
              ncpan=j+1
              ncpen=ncoup+1
              do jsh=ncoup,j+1,-1
                alpd(jsh+1)=alpd(jsh)
                jupc(jsh+1)=jupc(jsh)
                lupc(jsh+1)=lupc(jsh)
                mupc(jsh+1)=mupc(jsh)
              end do
              alpd(ncpan)=alpd(j)-zrad
              jupc(ncpan)=jupc(j)
              lupc(ncpan)=lupc(j)
              mupc(ncpan)=mupc(j)
C  DONE
CPB
c slmod end
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
C           IF (NLTRC) THEN
C             WRITE (6,*) ' TRACK TERMINATED'
C             WRITE (6,*) ' ITEST,ILIIN ',ITEST,ILIIN(IN)
C           ENDIF
            NCOUPE=J
            MPSURF=LUPC(NCOUPE)
            IPOLGN=LUPC(NCOUPE)
            NINCY=MUPC(NCOUPE)
            ZRAD=TSS
c slmod begin - debug - tr
            IF (printopt.GE.1.AND.printopt.LE.10)
     .        WRITE(6,'(5X,A,G12.5)')
     .          'TIMEP: *** TRACK ENDS ON NON-DEFAULT SURF *** (ZRAD)',
     .          zrad
c slmod end
            ISRFCL=0
            NINCX=0
            MRSURF=0
            MTSURF=0
            MASURF=0
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
c slmod begin - tr
c...note: Just trying to kill particle.  This problem is caused by inaccuracies
c         in the FindPoloidalCell routine.
          WRITE(0,*) 'TIMEP: ATTEMPTING TO KILL PARTICLE (AVOIDING'//
     .               ' LEARC1)'
          WRITE(6,*) 'TIMEP: ATTEMPTING TO KILL PARTICLE NPANU = ',NPANU
          ZRAD=1.0E+30
          RETURN
c
c          X0T=X00+VELX*ZRAD
c          Y0T=Y00+VELY*ZRAD
c          NN=LEARC1(X0T,Y0T,Z0T,IPOLGN,
c     .              NRCELL,NRCELL,.FALSE.,.FALSE.,
c     .              NPANU,'TIMEP       ')
c slmod end
        ENDIF
C
        MPSURF=0
        NINCY=0
C
311     CONTINUE
C       IF (NLTRC) THEN
C         WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C         WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C         WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C         WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C       ENDIF
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
c slmod begin - debug - tr
c            IF (printopt.GE.1.AND.printopt.LE.10) THEN
            IF (debugopt.NE.0) THEN
              WRITE(6,*)
              WRITE(0,'(5X,A)') 'TIMEP: ERROR: LEARC2 TIMEP CALL WILL'//
     .                          ' FAIL'
              WRITE(0,*) '   NPANU = ',npanu
            ENDIF
c slmod end
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
C       IF (NLTRC) THEN
C         WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C         WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C         WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C         WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C       ENDIF
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
          ALPD(NCOUP)=F-TWIN1
          TWIN1=F
C
          J1=J1+NINCY
          IF (J1.EQ.0) J1=NP2NDM
          IF (J1.EQ.NP2ND) J1=1
          IF (J1.NE.NJC) GOTO 100
C
C   LAST THETA-CELL
C
150       NCOUP=NCOUP+1
          JUPC(NCOUP)=NJC
          ALPD(NCOUP)=ZRAD-TWIN1
          X00=X000
          Y00=Y000
          NPCELL=NJC
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
c...sltmp
          STOP 'TIMER: NON-GENERALIZED GEOMETRY 01'
          DO 112 J=1,NP2ND
            V1=(YPOL(NRCELL,J)-Y00)*VELX-(XPOL(NRCELL,J)-X00)*VELY
            V2=(YPOL(NRCLLP,J)-Y00)*VELX-(XPOL(NRCLLP,J)-X00)*VELY
            LCUTY(J)=V1*V2.LE.0.
C           IF (NLTRC) THEN
C             IF (LCUTY(J)) WRITE (6,*) 'LCUTY=TRUE FOR ',J
C           ENDIF
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
c...sltmp
            STOP 'TIMER: UNGENERALIZED GEOMETRY 02'
            T1=((XPOL(NRCELL,I)-X00)*VVTY(NRCELL,I)-
     .          (YPOL(NRCELL,I)-Y00)*VVTX(NRCELL,I))
     .         /(VELX*VVTY(NRCELL,I)-VELY*VVTX(NRCELL,I)+EPS60)
C           IF (NLTRC) WRITE (6,*) 'I,T1 ',I,T1
            IF (T1.LT.0..OR.T1.GE.ZRAD) GOTO 114
C           IF (NLTRC) WRITE (6,*) 'VALID INTERSECTION AT T1= ',T1
            NCOUP=NCOUP+1
            LUPC(NCOUP)=I
            MUPC(NCOUP)=SIGN(1.D0,VELX*PPLNX(NRCELL,I)+
     .                            VELY*PPLNY(NRCELL,I))
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
C       IF (NLTRC.AND.NCOUP.GT.0) THEN
C         WRITE (6,*) ' NACH SORTIEREN '
C         WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C         WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C         WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C         WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C       ENDIF
C
        DO 125 J=1,NCOUP-1
          IF (ABS(ALPD(J+1)-ALPD(J)).LE.EPS30) THEN
            IF (JUPC(J).LE.0.OR.JUPC(J).GE.NP2ND) THEN
C             IF (NLTRC) THEN
C               WRITE (6,*) ' VERTAUSCHE ALPD(',J,') UND (',J+1,')'
C               WRITE (6,*) ' ALPD = ',ALPD(J),ALPD(J+1)
C               WRITE (6,*) ' JUPC = ',JUPC(J),JUPC(J+1)
C             ENDIF
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
          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
C
C  TRACK ENDS ON ONE OF THE NON DEFAULT POLOIDAL SURFACES
C
C           IF (NLTRC) THEN
C             WRITE (6,*) ' TRACK TERMINATED'
C             WRITE (6,*) ' ITEST,ILIIN ',ITEST,ILIIN(IN)
C           ENDIF
            NCOUPE=J
            MPSURF=LUPC(NCOUPE)
            IPOLGN=LUPC(NCOUPE)
            NINCY=MUPC(NCOUPE)
            ZRAD=TSS
            ISRFCL=0
            NINCX=0
            MRSURF=0
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
C         IF (NLTRC) THEN
C           WRITE (6,*) ' ALPD ',(ALPD(J),J=1,NCOUP)
C           WRITE (6,*) ' JUPC ',(JUPC(J),J=1,NCOUP)
C           WRITE (6,*) ' LUPC ',(LUPC(J),J=1,NCOUP)
C           WRITE (6,*) ' MUPC ',(MUPC(J),J=1,NCOUP)
C         ENDIF
C
c slmod begin - debug - tr
          IF (printopt.GE.1.AND.printopt.LE.10)
     .      WRITE(6,'(5X,A,2F10.4)')
     .        'TIMEP: Updating position 7 (VELX,VELY) ',VELX,VELY
c slmod end
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
C           IF (NLTRC) WRITE (6,*) 'MPTEST,F,DY ',MPTEST,F,DY
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
C         IF (NLTRC) THEN
C           WRITE (6,*) 'REENTRY FOUND, MPSURF,ZRAD = ',MPSURF,ZRAD
C           WRITE (6,*) 'IRCELL ',IRCELL
C         ENDIF
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
            Y0TEST=Y00+VELY*ZRD
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
C         IF (NLTRC) THEN
C           WRITE (6,*) 'NO REENTRY FOUND '
C         ENDIF
          Y00=Y00+ZRD*VELY
          GOTO 5000
        ENDIF
C
C  PARTICLE IN STANDARD MESH, RADIAL CELL NRCELL
C
2900    CONTINUE
        Y000=Y00+VELY*ZRAD
C       IF (NLTRC) WRITE (6,*) 'Y000,Y00,ZRAD ',Y000,Y00,ZRAD
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
          WRITE (6,*) 'ERROR IN TIMEP ',J2,J1,VELY
c slmod begin - tr
          WRITE (0,*) 'ERROR 1 IN TIMEP '
          ZRAD=1.0E+30
          RETURN
c
c          CALL EXIT
c slmod end
        ENDIF
C  TIME FROM Y00 TO PSURF
        IF (MPSURF.EQ.J2) THEN
          J1=J1+NINCY
          J2=J1+INCY
          GOTO 3000
        ENDIF
        DY=(PSURF(J2)-Y00)
        BB=VELY+EPS60
        F=DY/BB
C       IF (NLTRC) WRITE (6,*) 'J2,F,DY ',J2,F,DY
        IF (F.LE.ZRAD) THEN
          JUPC(NCOUP)=J1
          ALPD(NCOUP)=F-DUM
          DUM=F
C  STOP HISTORY AT NON DEFAULT STANDARD SURFACE J2
          ITEST=INMP2I(0,J2,0)
          IN=ITEST+NLIM
          IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
            ZRAD=F
            ISRFCL=0
            NPCELL=J1
            MPSURF=J2
            MRSURF=0
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
C       IF (NLTRC) WRITE (6,*) 'LAST CELL ',ZRAD,DUM
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
          WRITE (6,*) 'ERROR DETECTED IN TIMEP '
          WRITE (6,*) 'J,IADD+J,ALPD,JUPC ',J,IADD+J,ALPD(J),JUPC(J)
c slmod begin - debug - tr
          WRITE (6,*) 'ERROR DETECTED IN TIMEP (NP2ND) ',NP2ND
c          WRITE (0,*) 'ERROR DETECTED IN TIMEP '
          DO I = 1,NCOUP
            WRITE(6,*) J,NCOUP,JUPC(J),ALPD(J),NP2ND
          ENDDO
c slmod end
        ENDIF
C       IF (NLTRC) WRITE (6,*) 'TIMEP ',
C    .      J+IADD,CLPD(J+IADD),NUPC(J+IADD),NCOUNP(J+IADD)
5100  CONTINUE
C
      IF (ICOUT.LT.NCOUT.AND.MPSURF.EQ.0) THEN
        ICOUT=ICOUT+1
        IPOLGN=0
        IF (ICOUT.EQ.NCOUT) IPOLGN=IPOLGS
        ZRAD=BLPD(ICOUT)
c slmod begin - tr (already there)
CPB
c... This is a bug fix sent by Detlev on March 16, 2001.
C  RESTORE POLOIDAL SECTIONS FOR THE REST OF THE TRACK
C  THIS IS NEEDED ONLY IN LEVGEO=3 OPTION. OTHERWISE NCPAN=NCPEN=0 
        jn=0
        do j=ncpan,ncpen
          jn=jn+1
          alpd(jn)=alpd(j)
          jupc(jn)=jupc(j)
          lupc(jn)=lupc(j)
          mupc(jn)=mupc(j)
        end do
	  ncoup=jn
C  DONE
CPB
c slmod end
        NTCELL=NCOUNT(ICOUT)
        IADD=IADD+NCOUP
C       IF (NLTRC)
C    .  WRITE (6,*) 'NEXT TOR. CELL: ZRAD,NTCELL,IADD ',ZRAD,NTCELL,IADD
        GOTO 10000
      ENDIF
C
      NCOU=IADD+NCOUP
C
      SUM=0.
      DO 5110 ICOU=1,NCOU
        SUM=SUM+CLPD(ICOU)
        NCOUNT(ICOU)=NCOUNS(ICOU)
C       WRITE (6,*) 'ICOU,NCOUNT,CLPD ',ICOU,NCOUNT(ICOU),CLPD(ICOU)
5110  CONTINUE
      IF (MPSURF.EQ.0.AND.ABS(SUM-ZRADS).GT.EPS10) THEN
        WRITE (6,*) 'ERROR IN TIMEP: NPANU,SUM,ZRADS ',
     .                               NPANU,SUM,ZRADS
c slmod begin - tr (already there)
        WRITE (6,*) 'ERROR 2 IN TIMEP '
        WRITE (0,*) 'ERROR 2 IN TIMEP '
        ZRAD=1.0E+30
        RETURN
c
c       CALL EXIT
c slmod end
      ENDIF
C
      ZRAD=SUM
c slmod begin - debug - tr
      IF (printopt.GE.1.AND.printopt.LE.10) THEN
        WRITE(6,'(5X,A,5I4,3X,1P,E12.5)')
     .    'TIMEP: Output (MPSURF NP,IP,IR,IPOLGN ZRAD) ',
     .    MPSURF,NPCELL,IPCELL,IRCELL,IPOLGN,ZRAD
        WRITE(6,'(5X,A)') 'TIMEP: '
      ENDIF
c slmod end
      RETURN
C
991   CONTINUE
      WRITE (6,*) 'REENTRANCE FROM VACUUM REGION AND LEVGEO=1 IN TIMEP'
      CALL EXIT
      END
C
      SUBROUTINE TIMET (ZRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
C    TO BE WRITTEN
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
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'CUPD'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
      INCLUDE 'CLOGAU'
      INCLUDE 'CLGIN'
      INCLUDE 'COMSOU'
C
C     CALCULATE TIME SEGMENTS FOR 2D-PROFILES,
C     RADIALLY AND TOROIDALLY RESOLVED
C
      ZRADS=ZRAD
C
      IF (.NOT.NLTRZ) GOTO 1000
C
C  CYLINDRICAL OR CARTHESIAN CO-ORDINATE SYSTEM
C
C     IF (NLTRC) THEN
C       WRITE (6,*) 'TIMET: ZRAD,NRCELL,NTCELL,MTSURF '
C       WRITE (6,*) '      ',ZRAD,NRCELL,NTCELL,MTSURF
C       WRITE (6,*) 'INITIAL: X0,Z01 ',X0,Z01
C     ENDIF
C
C  PARTICLE OUTSIDE STANDARD MESH ?
C
      IF (NRCELL.EQ.0) THEN
C
C
C  PARTICLE OUTSIDE STANDARD MESH
C  CHECK AT NON DEFAULT POLOIDAL SURFACES
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
C           IF (NLTRC) WRITE (6,*) 'MTTEST,F,DZ ',MTTEST,F,DZ
            IF (F.LE.ZRD.AND.F.GT.0.D0) THEN
              X0TEST=X0+VELX*F
              IF (X0TEST.GE.RSURF(1).AND.X0TEST.LE.RSURF(NR1ST)) THEN
                IRSAVE=LEARCA(X0TEST,RSURF,1,NR1ST,1,'TIMET 1    ')
                NCOUT=1
                JUPC(NCOUT)=1
                MTSAVE=MTTEST
                ZRD=F
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
C         IF (NLTRC) THEN
C           WRITE (6,*) 'REENTRY FOUND, MTSURF,ZRAD = ',MTSURF,ZRAD
C           WRITE (6,*) 'IRCELL ',IRCELL
C         ENDIF
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
            Z0TEST=Z00+VELZ*ZRD
            IF (Z0TEST.GE.ZSURF(1).AND.Z0TEST.LE.ZSURF(NT3RD)) THEN
              ITCELL=LEARCA(Z0TEST,ZSURF,1,NT3RD,1,'TIMET 2     ')
            ELSE
              MRSURF=0
              NINCX=0
            ENDIF
          ENDIF
          MTSURF=0
          NINCZ=0
C         IF (NLTRC) THEN
C           WRITE (6,*) 'NO REENTRY FOUND '
C         ENDIF
          Z01=Z01+ZRD*VELZ
          GOTO 5000
        ENDIF
C
      ENDIF
C
C  PARTICLE IN STANDARD MESH, RADIAL CELL NRCELL
C
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
      BB=VELZ+EPS60
      F=DZ/BB
      IF (F.LE.ZRAD) THEN
        KUPC(NCOUT)=J1
        BLPD(NCOUT)=F-DUM
        DUM=F
C  STOP HISTORY AT NON DEFAULT STANDARD SURFACE J2
        ITEST=INMP3I(0,0,J2)
        IN=ITEST+NLIM
        IF (ITEST.NE.0.AND.ILIIN(IN).NE.0) THEN
          ZRAD=F
          ISRFCL=0
          NTCELL=J1
          MTSURF=J2
          MRSURF=0
          MASURF=0
          IPOLGN=0
          ITCELL=NTCELL
          Z01=ZSURF(J2)
C         IF (NLTRC) WRITE (6,*) 'STOP AT MTSURF ',MTSURF
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
      KUPC(NCOUT)=J1
      BLPD(NCOUT)=ZRAD-DUM
      ITCELL=NTCELL
      Z01=Z001
C
      GOTO 5000
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN TIMET, Z SURFACE INDEX OUT OF RANGE  '
      WRITE (6,*) 'NPANU,Z0,Z01,ZRAD,VELZ,NTCELL '
      WRITE (6,*)  NPANU,Z0,Z01,ZRAD,VELZ,NTCELL
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
C     IF (NLTRC) THEN
C       WRITE (6,*) 'TIMET: ZRAD,NRCELL,NTCELL,MTSURF '
C       WRITE (6,*) '      ',ZRAD,NRCELL,NTCELL,MTSURF
C       WRITE (6,*) 'INITIAL: X01,Z01,PHI ',X01,Z01,PHI
C     ENDIF
C
      NERR=0
      NCOUT=1
      KUPC(1)=NTCELL
C
C     IF (NLSRFZ) THEN ....
      NLSRFZ=.FALSE.
c slmod begin - tr
c...  This may not be required any more.  The problem this addition addressed
c     was related to always launching the target neutral at the segment boundary,
c     which is not done any more.
      NLSRFT=.FALSE.
c slmod end
C
1010  CONTINUE
C  PHI0 IS THE PHI AT THE CENTER OF THE CURRENT TOROIDAL CELL
      PHI0=PHI-ATAN2(Z01,X01)
C
C     IF (NLTRC) THEN
C       TTT=Z01/(X01*TANAL)
C       IF (ABS(TTT).GT.1.+EPS10) THEN
C         WRITE (6,*) 'NPANU ',NPANU
C         WRITE (6,*) 'X01,Z01 OUT OF RANGE IN TIMET'
C         WRITE (6,*) X01,Z01,TTT
C         CALL EXIT
C       ENDIF
C     ENDIF
C
      Z001=Z01+ZRAD*VELZ
      X001=X01+ZRAD*VELX
      IF (ZRAD.LT.1.D30.AND.X01*X001.GT.0.D0) THEN
        ITT=IDINT(Z001/(X001*TANAL))
C       IF (NLTRC) WRITE (6,*) 'TIMET 1 ',X01,Z01,X001,Z001,ITT
      ELSE
        TO=(Z01-X01*TANAL)/(TANAL*VELX-VELZ)
        XTO=X01+TO*VELX
        TU=(Z01+X01*TANAL)/(-TANAL*VELX-VELZ)
        XTU=X01+TU*VELX
C       IF (NLTRC) WRITE (6,*) 'TU,TO ',TU,TO,XTU,XTO
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
C       IF (NLTRC) WRITE (6,*) 'TIMET 2 ',X01,Z01,ITT
      ENDIF
C
c slmod begin - tr
      NINCS=0
c slmod end
      IF (ITT) 1100,1200,1300
C
C  NO INTERSECTION WITH TOROIDAL SURFACE
C
1200  CONTINUE
      NNTCLL=0
      MTSURF=0
      Z01=Z001
      X01=X001
C  IN CASE ZRAD=1.D30, IS NEXT STATEMENT IS NONSENSE, BUT CORRECTED
C                      FOR IN SUBR. STDCOL, WHICH MUST BE CALLED NEXT
C                      FOR A POLOIDAL SURFACE (OTHERWISE: ERROR EXIT)
      PHI=PHI0+ATAN2(Z01,X01)
      BLPD(1)=ZRAD
C     IF (NLTRC) WRITE (6,*) 'FINAL 1: X01,Z01,PHI ',X01,Z01,PHI
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
C         IF (NLTRC) WRITE (6,*) 'TRY AGAIN IN TIMET'
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
c slmod begin - tr
      NINCS=1
c slmod end
      IF (NNTCLL.GE.NTTRA) NNTCLL=1
      MTSURF=NNTCLL
C     IF (NLTRC) THEN
C       WRITE (6,*) 'FINAL 2: X01,Z01,PHI ',X01,Z01,PHI
C       WRITE (6,*) 'ZRAD,ISRFCL ',ZRAD,ISRFCL
C     ENDIF
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
C         IF (NLTRC) WRITE (6,*) 'TRY AGAIN IN TIMET'
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
c slmod begin - tr
      NINCS=-1
c slmod end
      IF (NNTCLL.LE.0) NNTCLL=NTTRAM
      MTSURF=NNTCLL+1
C     IF (NLTRC) THEN
C       WRITE (6,*) 'FINAL 3: X01,Z01,PHI ',X01,Z01,PHI
C       WRITE (6,*) 'ZRAD,ISRFCL ',ZRAD,ISRFCL
C     ENDIF
      GOTO 5000
C
C  DISCRETE TOROIDAL APPROXIMATION FINISHED
C
C
5000  CONTINUE
C     IF (NLTRC) WRITE (6,*) 'NCOUT= ',NCOUT
      DO 5100 J=1,NCOUT
        CLPD(J)=BLPD(J)
        NUPC(J)=(KUPC(J)-1)*NP2T3
        NCOUNT(J)=KUPC(J)
        IF (CLPD(J).LE.0..OR.KUPC(J).LE.0.OR.
     .      (KUPC(J).GE.NT3RD.AND.NLTOR)) THEN
          WRITE (6,*) 'ERROR DETECTED IN TIMET '
          WRITE (6,*) 'J,BLPD,KUPC ',J,BLPD(J),KUPC(J)
        ENDIF
C       IF (NLTRC) THEN
C         WRITE (6,*) 'TIMET: J,BLPD,NUPC,NCOUNT ',
C    .                        J,BLPD(J),NUPC(J),NCOUNT(J)
C       ENDIF
5100  CONTINUE
C     IF (NLTRC) THEN
C       WRITE (6,*) 'MTSURF,NLSRFZ,NINCZ,IRCELL,IPCELL '
C       WRITE (6,*)  MTSURF,NLSRFZ,NINCZ,IRCELL,IPCELL
C       WRITE (6,*) 'INMP3I ',INMP3I(IRCELL,IPCELL,MTSURF)
C     ENDIF
C
      NCOU=NCOUT
C
      SUM=0.
      DO 5110 ICOU=1,NCOU
        SUM=SUM+CLPD(ICOU)
C       WRITE (6,*) 'ICOU,NCOUNT,CLPD ',ICOU,NCOUNT(ICOU),CLPD(ICOU)
5110  CONTINUE
      IF (MTSURF.EQ.0.AND.ABS(SUM-ZRADS).GT.EPS10) THEN
        WRITE (6,*) 'ERROR IN TIMET: NPANU,SUM,ZRADS ',
     .                               NPANU,SUM,ZRADS
      ENDIF
      ZRAD=SUM
C
      RETURN
C
9998  CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN TIMET. RETURN ZRAD=1.D30'
      WRITE (6,*) 'NPANU,AA,BB ',NPANU,AA,BB
      RETURN
9999  CONTINUE
      WRITE (6,*) 'INVALID OPTION IN TIMET. EXIT CALLED '
      CALL EXIT
      END
C
      SUBROUTINE XSHADD (XSH,ILINI,ILEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  SHIFT CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES IN X DIRECTION
C  OLD ORIGIN: X0=0.    (X,Y,Z)  SYSTEM
C  NEW ORIGIN: XN=-XSH   (X',Y,Z) SYSTEM
C
C     X'=X+XSH
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'
C
      XQ=XSH*XSH
C
      DO 100 I=ILINI,ILEND
C
C
C   CHANGE COEFFICIENTS FOR ALGEBRAIC EQUATION
          A0LM(I)=A0LM(I)-A1LM(I)*XSH+A4LM(I)*XQ
          A1LM(I)=A1LM(I)-2.*A4LM(I)*XSH
          A2LM(I)=A2LM(I)-A7LM(I)*XSH
          A3LM(I)=A3LM(I)-A8LM(I)*XSH
          IF (RLB(I).LT.0.) THEN
C   CHANGE COEFFICIENTS IN LINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 10 J=1,ILIN(I)
              ALIMS(I,J)=ALIMS(I,J)-XLIMS(I,J)*XSH
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
              ALIMS0(I,J)=ALIMS0(I,J)-XLIMS1(I,J)*XSH+XLIMS2(I,J)*XQ
              XLIMS1(I,J)=XLIMS1(I,J)-2.*XLIMS2(I,J)*XSH
              YLIMS1(I,J)=YLIMS1(I,J)-XLIMS3(I,J)*XSH
              ZLIMS1(I,J)=ZLIMS1(I,J)-YLIMS3(I,J)*XSH
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER
            XLIMS1(I,1)=XLIMS1(I,1)+XSH
            XLIMS2(I,1)=XLIMS2(I,1)+XSH
C
C
          ELSE
C   POINT OPTIONS
            P1(1,I)=P1(1,I)+XSH
            P2(1,I)=P2(1,I)+XSH
            IF (RLB(I).GE.3.) P3(1,I)=P3(1,I)+XSH
            IF (RLB(I).GE.4.) P4(1,I)=P4(1,I)+XSH
            IF (RLB(I).GE.5.) P5(1,I)=P5(1,I)+XSH
            IF (RLB(I).GE.6.) P6(1,I)=P6(1,I)+XSH
          ENDIF
C
100   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE ZSHADD (ZSH,ILINI,ILEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  SHIFT CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES IN Z DIRECTION
C  OLD ORIGIN: ZO=O.    (X,Y,Z)  SYSTEM
C  NEW ORIGIN: ZN=-ZSH   (X,Y,Z') SYSTEM
C
C    Z'=Z+ZSH
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'
C
      ZQ=ZSH*ZSH
C
      DO 100 I=ILINI,ILEND
C
C
C   CHANGE COEFFICIENTS FOR ALGEBRAIC EQUATION
          A0LM(I)=A0LM(I)-A3LM(I)*ZSH+A6LM(I)*ZQ
          A3LM(I)=A3LM(I)-2.*A6LM(I)*ZSH
          A1LM(I)=A1LM(I)-A8LM(I)*ZSH
          A2LM(I)=A2LM(I)-A9LM(I)*ZSH
          IF (RLB(I).LT.0.) THEN
C   CHANGE COEFFICIENTS IN LINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 10 J=1,ILIN(I)
              ALIMS(I,J)=ALIMS(I,J)-ZLIMS(I,J)*ZSH
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
              ALIMS0(I,J)=ALIMS0(I,J)-ZLIMS1(I,J)*ZSH+ZLIMS2(I,J)*ZQ
              ZLIMS1(I,J)=ZLIMS1(I,J)-2.*ZLIMS2(I,J)*ZSH
              XLIMS1(I,J)=XLIMS1(I,J)-YLIMS3(I,J)*ZSH
              YLIMS1(I,J)=YLIMS1(I,J)-ZLIMS3(I,J)*ZSH
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER
            ZLIMS1(I,1)=ZLIMS1(I,1)+ZSH
            ZLIMS2(I,1)=ZLIMS2(I,1)+ZSH
C
          ELSE
C
C   POINT OPTIONS
            P1(3,I)=P1(3,I)+ZSH
            P2(3,I)=P2(3,I)+ZSH
            IF (RLB(I).GE.3.) P3(3,I)=P3(3,I)+ZSH
            IF (RLB(I).GE.4.) P4(3,I)=P4(3,I)+ZSH
            IF (RLB(I).GE.5.) P5(3,I)=P5(3,I)+ZSH
            IF (RLB(I).GE.6.) P6(3,I)=P6(3,I)+ZSH
          ENDIF
C
100   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE YSHADD (YSH,ILINI,ILEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  SHIFT CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES IN Y DIRECTION
C  OLD ORIGIN: YO=O.    (X,Y,Z)  SYSTEM
C  NEW ORIGIN: YN=-YSH   (X,Y',Z) SYSTEM
C
C    Y'=Y+YSH
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'
C
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
              ALIMS(I,J)=ALIMS(I,J)-YLIMS(I,J)*YSH
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
              ALIMS0(I,J)=ALIMS0(I,J)-YLIMS1(I,J)*YSH+YLIMS2(I,J)*YQ
              YLIMS1(I,J)=YLIMS1(I,J)-2.*YLIMS2(I,J)*YSH
              XLIMS1(I,J)=XLIMS1(I,J)-XLIMS3(I,J)*YSH
              ZLIMS1(I,J)=ZLIMS1(I,J)-ZLIMS3(I,J)*YSH
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER
            YLIMS1(I,1)=YLIMS1(I,1)+YSH
            YLIMS2(I,1)=YLIMS2(I,1)+YSH
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
C
C
      SUBROUTINE ROTADD(A,AI,ILINI,ILEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  TRANSFORM CO-ORDINATE SYSTEM FOR ADDITIONAL SURFACES
C  CO-ORDINATES OF A POINT IN OLD SYSTEM (X,Y,Z)
C  CO-ORDINATES OF A POINT IN NEW SYSTEM (X',Y',Z')
C
C    X'         X     X          X'                           -1  T
C    Y' =   A * Y  ;  Y  =  AI * Y' ; SUBROUTINE ASSUMES: AI=A  =A
C    Z'         Z     Z          Z'
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CADGEO'
      INCLUDE 'CTRCEI'
      DIMENSION A(3,3),AI(3,3)
C
C
      DO 100 I=ILINI,ILEND
C
C   CHANGE COEFFICIENTS FOR ALGEBRAIC EQUATION, TRY TO KEEP INVARIANTS
          S=A4LM(I)+A5LM(I)+A6LM(I)
          DEL=DETER(A4LM(I),A7LM(I)/2.,A9LM(I)/2.,A7LM(I)/2.,
     .              A5LM(I),A8LM(I)/2.,A9LM(I)/2.,A8LM(I)/2.,
     .              A6LM(I))
          T=A5LM(I)*A6LM(I)+A6LM(I)*A4LM(I)+A4LM(I)*A5LM(I)-
     .     (A9LM(I)**2+A8LM(I)**2+A7LM(I)**2)*0.25
C         IF (TRCPLT) THEN
C           WRITE (6,*) 'INVARIANTS FROM ROTADD: I= ',I
C           WRITE (6,*) 'BEFORE: S,DEL,T= ',S,DEL,T
C           WRITE (6,*) ' A0...A9 ',A0LM(I),A1LM(I),A2LM(I),A3LM(I),
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
          DEL=DETER(A4LM(I),A7LM(I)/2.,A9LM(I)/2.,A7LM(I)/2.,
     .              A5LM(I),A8LM(I)/2.,A9LM(I)/2.,A8LM(I)/2.,
     .              A6LM(I))
          T=A5LM(I)*A6LM(I)+A6LM(I)*A4LM(I)+A4LM(I)*A5LM(I)-
     .     (A9LM(I)**2+A8LM(I)**2+A7LM(I)**2)*0.25
C         WRITE (6,*) 'AFTER:  S,DEL,T= ',S,DEL,T
C
1         IF (RLB(I).LT.0.) THEN
C   CHANGE COEFFICIENTS IN LINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 10 J=1,ILIN(I)
C             ALIMS(I,J)=ALIMS(I,J)
              A1S=XLIMS(I,J)*AI(1,1)+YLIMS(I,J)*AI(2,1)+
     .            ZLIMS(I,J)*AI(3,1)
              A2S=XLIMS(I,J)*AI(1,2)+YLIMS(I,J)*AI(2,2)+
     .            ZLIMS(I,J)*AI(3,2)
              A3S=XLIMS(I,J)*AI(1,3)+YLIMS(I,J)*AI(2,3)+
     .            ZLIMS(I,J)*AI(3,3)
              XLIMS(I,J)=A1S
              YLIMS(I,J)=A2S
              ZLIMS(I,J)=A3S
10          CONTINUE
C   CHANGE COEFFICIENTS IN NONLINEAR INEQUALITIES BOUNDING THE SURFACE
            DO 20 J=1,ISCN(I)
C             ALIMS0(I,J)=ALIMS0(I,J)
              A1S=XLIMS1(I,J)*AI(1,1)+YLIMS1(I,J)*AI(2,1)+
     .            ZLIMS1(I,J)*AI(3,1)
              A2S=XLIMS1(I,J)*AI(1,2)+YLIMS1(I,J)*AI(2,2)+
     .            ZLIMS1(I,J)*AI(3,2)
              A3S=XLIMS1(I,J)*AI(1,3)+YLIMS1(I,J)*AI(2,3)+
     .            ZLIMS1(I,J)*AI(3,3)
              XLIMS1(I,J)=A1S
              YLIMS1(I,J)=A2S
              ZLIMS1(I,J)=A3S
C
              S=XLIMS2(I,J)+YLIMS2(I,J)+ZLIMS2(I,J)
              A4S=XLIMS2(I,J)*AI(1,1)**2+YLIMS2(I,J)*AI(2,1)**2+
     .            ZLIMS2(I,J)*AI(3,1)**2+
     .            XLIMS3(I,J)*AI(1,1)*AI(2,1)+
     .            YLIMS3(I,J)*AI(1,1)*AI(3,1)+
     .            ZLIMS3(I,J)*AI(2,1)*AI(3,1)
              A5S=XLIMS2(I,J)*AI(1,2)**2+YLIMS2(I,J)*AI(2,2)**2+
     .            ZLIMS2(I,J)*AI(3,2)**2+
     .            XLIMS3(I,J)*AI(1,2)*AI(2,2)+
     .            YLIMS3(I,J)*AI(1,2)*AI(3,2)+
     .            ZLIMS3(I,J)*AI(2,2)*AI(3,2)
              A6S=S-A4S-A5S
C
              A7S=XLIMS2(I,J)*(AI(1,1)*AI(1,2)+AI(1,1)*AI(1,2))+
     .            YLIMS2(I,J)*(AI(2,1)*AI(2,2)+AI(2,1)*AI(2,2))+
     .            ZLIMS2(I,J)*(AI(3,1)*AI(3,2)+AI(3,1)*AI(3,2))+
     .            XLIMS3(I,J)*(AI(1,1)*AI(2,2)+AI(1,2)*AI(2,1))+
     .            YLIMS3(I,J)*(AI(1,1)*AI(3,2)+AI(1,2)*AI(3,1))+
     .            ZLIMS3(I,J)*(AI(2,1)*AI(3,2)+AI(2,2)*AI(3,1))
              A8S=XLIMS2(I,J)*(AI(1,1)*AI(1,3)+AI(1,1)*AI(1,3))+
     .            YLIMS2(I,J)*(AI(2,1)*AI(2,3)+AI(2,1)*AI(2,3))+
     .            ZLIMS2(I,J)*(AI(3,1)*AI(3,3)+AI(3,1)*AI(3,3))+
     .            XLIMS3(I,J)*(AI(1,1)*AI(2,3)+AI(1,3)*AI(2,1))+
     .            YLIMS3(I,J)*(AI(1,1)*AI(3,3)+AI(1,3)*AI(3,1))+
     .            ZLIMS3(I,J)*(AI(2,1)*AI(3,3)+AI(2,3)*AI(3,1))
              A9S=XLIMS2(I,J)*(AI(1,2)*AI(1,3)+AI(1,2)*AI(1,3))+
     .            YLIMS2(I,J)*(AI(2,2)*AI(2,3)+AI(2,2)*AI(2,3))+
     .            ZLIMS2(I,J)*(AI(3,2)*AI(3,3)+AI(3,2)*AI(3,3))+
     .            XLIMS3(I,J)*(AI(2,1)*AI(2,3)+AI(1,3)*AI(2,2))+
     .            YLIMS3(I,J)*(AI(1,2)*AI(3,3)+AI(1,3)*AI(3,2))+
     .            ZLIMS3(I,J)*(AI(2,2)*AI(3,3)+AI(2,3)*AI(3,2))
              XLIMS2(I,J)=A4S
              YLIMS2(I,J)=A5S
              ZLIMS2(I,J)=A6S
              XLIMS3(I,J)=A7S
              YLIMS3(I,J)=A8S
              ZLIMS3(I,J)=A9S
20          CONTINUE
          ELSEIF (RLB(I).GT.0.AND.RLB(I).LT.2) THEN
C  BOUNDED BY QUADER: CHANGE TO RLB(I)=-6 OPTION AND DEFINE 6 LINEAR
C                     INEQUALITIES
            ALIMS(I,1)=XLIMS1(I,1)
            XLIMS(I,1)=-1.
            YLIMS(I,1)=0.
            ZLIMS(I,1)=0.
            ALIMS(I,2)=-XLIMS2(I,1)
            XLIMS(I,2)=1.
            YLIMS(I,2)=0.
            ZLIMS(I,2)=0.
            ALIMS(I,3)=YLIMS1(I,1)
            XLIMS(I,3)=0.
            YLIMS(I,3)=-1.
            ZLIMS(I,3)=0.
            ALIMS(I,4)=-YLIMS2(I,1)
            XLIMS(I,4)=0.
            YLIMS(I,4)=1.
            ZLIMS(I,4)=0.
            ALIMS(I,5)=ZLIMS1(I,1)
            XLIMS(I,5)=0.
            YLIMS(I,5)=0.
            ZLIMS(I,5)=-1.
            ALIMS(I,6)=-ZLIMS2(I,1)
            XLIMS(I,6)=0.
            YLIMS(I,6)=0.
            ZLIMS(I,6)=1.
            IF (RLB(I).EQ.1.5) THEN
              WRITE (6,*) 'ROTADD: TO BE WRITTEN, RLB=1.5'
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
C
C
      SUBROUTINE SNEIGH
C  DETERMINE THE INDICES OF THE FOUR NEIGHBORING CELLS OF
C  AN EIRENE CELL    : NGHPOL
C  AN EIRENE SURFCASE: NGHPLS
C  SIDE NUMBERING:
C
C       (IR+1,IP+1)          (IR,IP+1)
C                       2
C                 +-----------+
C                 |           |
C                 |           |
C               3 |           | 1
C                 |           |
C                 |           |
C                 +-----------+
C                       4
C         (IR+1,IP)          (IR,IP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'CGEOM'
      INCLUDE 'CCONA'
      INCLUDE 'CPOLYG'

      DO IR=1,NR1ST
        DO IP=1,NP2ND
          NGHPOL(1,IR,IP)=0
          NGHPOL(2,IR,IP)=0
          NGHPOL(3,IR,IP)=0
          NGHPOL(4,IR,IP)=0
          NGHPLS(2,IR,IP)=0
          NGHPLS(4,IR,IP)=0
        ENDDO
      ENDDO

      DO IR=1,NR1STM
        DO K=1,NPPLG
          DO IP=NPOINT(1,K),NPOINT(2,K)
C  RADIAL NEIGHBORS
c slmod begin - new
c...        Not sure why NGHPOL(3,... used the remainder function.  It is
c           not referred to anywhere in this version, but perhaps in 
c           EIRENE02 it is more of an issue:
            NGHPOL(1,IR,IP)=IR-1
            NGHPOL(3,IR,IP)=IR+1

c            WRITE(0,'(A,4I6)') 'POL:',ir,ip,
c     .                 NGHPOL(3,IR,IP),NGHPOL2(2,IR,IP)

            NGHPOL(1,IR,IP)=NGHPOL2(1,IR,IP)
            NGHPOL(3,IR,IP)=NGHPOL2(2,IR,IP)
 
c
c            NGHPOL(1,IR,IP)=IR-1
c            NGHPOL(3,IR,IP)=MOD(IR+1,NR1ST)
c slmod end
C  POLOIDAL NEIGHBORS
            IF (IP.GT.NPOINT(1,K).AND.IP.LT.NPOINT(2,K)-1) THEN
C  INNER CELL
              NGHPOL(4,IR,IP)=IP-1
              NGHPOL(2,IR,IP)=IP+1
              NGHPLS(4,IR,IP)=IP-1
              NGHPLS(2,IR,IP)=IP
            ELSE IF (IP.EQ.NPOINT(1,K)) THEN
C  FIRST CELL OR SURFACE IN A PART OF A POLYGON
              NGHPOL(2,IR,IP)=IP+1
              NGHPLS(2,IR,IP)=IP
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
              DO IPART=1,NPPLG
                JP=NPOINT(2,IPART)
c slmod begin - tr - new
c DO NOT LINK THE FIRST POLOIDAL SURFACE, AND ALSO STOP
c CELLS LINKING TO ANOTHER CELL IN THE SAME GROUP
                IF (K.EQ.1.OR.
     .              (K.EQ.NPPLG.AND.K.EQ.IPART)) CYCLE
                IF (GRIDOPT.EQ.1) THEN
                  EP=EPS10
                  IF (((XVERT(IR,IP,1)-XVERT(IR,JP,1))**2+
     .                 (YVERT(IR,IP,1)-YVERT(IR,JP,1))**2).LT.EP.AND.
     .                ((XVERT(IR,IP,2)-XVERT(IR,JP,2))**2+
     .                 (YVERT(IR,IP,2)-YVERT(IR,JP,2))**2).LT.EP) THEN
                    NGHPOL(4,IR,IP)=JP-1
                    NGHPLS(4,IR,IP)=JP-1
                    GOTO 1
                  ENDIF
                ELSE
                  IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
     .                 (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS10 .AND.
     .                ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .                 (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS10) THEN
                    NGHPOL(4,IR,IP)=JP-1
                    NGHPLS(4,IR,IP)=JP-1
                    GOTO 1
                  ENDIF
                ENDIF
c
c                IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
c     .               (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS10 .AND.
c     .              ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
c     .               (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS10) THEN
c                  NGHPOL(4,IR,IP)=JP-1
c                  NGHPLS(4,IR,IP)=JP-1
c                  GOTO 1
c                ENDIF
c slmod end
              ENDDO
            ELSEIF (IP.EQ.NPOINT(2,K)-1) THEN
C  LAST CELL IN A PART OF A POLYGON
              NGHPOL(4,IR,IP)=IP-1
              NGHPLS(4,IR,IP)=IP-1
              NGHPLS(2,IR,IP)=IP
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
              DO IPART=1,NPPLG
                JP=NPOINT(1,IPART)
c slmod begin - tr - new
                IF ( K.EQ.NPPLG.OR.
     .              (K.EQ.1.AND.K.EQ.IPART)) CYCLE
                IF (GRIDOPT.EQ.1) THEN
                  EP=EPS10
                  IF (((XVERT(IR,IP+1,1)-XVERT(IR,JP,1))**2+
     .                 (YVERT(IR,IP+1,1)-YVERT(IR,JP,1))**2).LT.EP.AND.
     .                ((XVERT(IR,IP+1,2)-XVERT(IR,JP,2))**2+
     .                 (YVERT(IR,IP+1,2)-YVERT(IR,JP,2))**2).LT.EP) THEN
                    NGHPOL(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ELSE
                  IF (((XPOL(IR,IP+1)-XPOL(IR,JP))**2+
     .                (YPOL(IR,IP+1)-YPOL(IR,JP))**2).LT.EPS10 .AND.
     .               ((XPOL(IR+1,IP+1)-XPOL(IR+1,JP))**2+
     .                (YPOL(IR+1,IP+1)-YPOL(IR+1,JP))**2).LT.EPS10) THEN
                    NGHPOL(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ENDIF
c
c                IF (((XPOL(IR,IP+1)-XPOL(IR,JP))**2+
c     .               (YPOL(IR,IP+1)-YPOL(IR,JP))**2).LT.EPS10 .AND.
c     .              ((XPOL(IR+1,IP+1)-XPOL(IR+1,JP))**2+
c     .               (YPOL(IR+1,IP+1)-YPOL(IR+1,JP))**2).LT.EPS10) THEN
c                  NGHPOL(2,IR,IP)=JP
c                  GOTO 1
c                ENDIF
c slmod end
              ENDDO
            ELSEIF (IP.EQ.NPOINT(2,K)) THEN
C  LAST SURFACE IN A PART OF A POLYGON
              NGHPLS(4,IR,IP)=IP-1
C  LOOK FOR EQUALITY WITH OTHER POLOIDAL POLYGONS
              DO IPART=1,NPPLG
                JP=NPOINT(1,IPART)
c slmod begin - new - tr - new
                IF (K.EQ.IPART) CYCLE
                IF (GRIDOPT.EQ.1) THEN
                  EP=EPS10
                  IF (((XVERT(IR,IP,1)-XVERT(IR,JP,1))**2+
     .                 (YVERT(IR,IP,1)-YVERT(IR,JP,1))**2).LT.EP.AND.
     .                ((XVERT(IR,IP,2)-XVERT(IR,JP,2))**2+
     .                 (YVERT(IR,IP,2)-YVERT(IR,JP,2))**2).LT.EP) THEN
                    NGHPLS(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ELSE
                  IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
     .                 (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS10 .AND.
     .                ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
     .                 (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS10) THEN
                    NGHPLS(2,IR,IP)=JP
                    GOTO 1
                  ENDIF
                ENDIF
c
c                IF (((XPOL(IR,IP)-XPOL(IR,JP))**2+
c     .               (YPOL(IR,IP)-YPOL(IR,JP))**2).LT.EPS10 .AND.
c     .              ((XPOL(IR+1,IP)-XPOL(IR+1,JP))**2+
c     .               (YPOL(IR+1,IP)-YPOL(IR+1,JP))**2).LT.EPS10) THEN
c                  NGHPLS(2,IR,IP)=JP
c                  GOTO 1
c                ENDIF
c slmod end
              ENDDO
            ENDIF
1         ENDDO
        ENDDO
      ENDDO
c slmod begin - not tr
      WRITE (6,*) 'CONNECTION MAP:'
      WRITE (6,*) '  IR    IP    S1    S2    S3    S4'
      DO IR=1,NR1STM
        DO IP=1,NP2NDM
          WRITE (6,'(8I6)') IR,IP,(NGHPOL(K,IR,IP),K=1,4),
     .             nghpol2(1,ir,ip),nghpol2(2,ir,ip)
c          WRITE (6,'(6I6)') IR,IP,(NGHPLS(K,IR,IP),K=1,4)
        ENDDO
      ENDDO
c slmod end
      END
