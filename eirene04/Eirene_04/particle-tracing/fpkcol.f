C
C
      SUBROUTINE FPKCOL(*,*)
C
C  FOKKER PLANCK ELASTIC COLLISION
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CFPLK
      USE CLOGAU
      USE CUPD
      USE CGRID
      USE CGEOM
      USE CZT1
      USE COMPRT
      USE CLGIN
      USE COUTAU
      USE COMXS

      IMPLICIT NONE

      REAL(DP) :: DUR, ENEW, VNEW, WS
      INTEGER :: IOLD, LEARC2, NCELLT
C  RETURN 1:  NOT IN USE
C  RETURN 2:  START COMPLETELY NEW TEST ION TRACK, SAME SPECIES
C
C  SAVE INCIDENT SPECIES: IOLD
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,7)
      IOLD=IION
      NCELLT=NCLTAL(NCELL)
C
      X0=X0+VLXPAR*ZT
      Y0=Y0+VLYPAR*ZT
      Z0=Z0+VLZPAR*ZT
      TIME=TIME+ZT/VELPAR
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=LEARC2(X0,Y0,NRCELL,NPANU,'FOLION 2     ')
      ELSEIF (NLFEM) THEN
        IPOLG=0
      ENDIF
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      MSURF=0
      IF (NLTRA) PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
C
255   CONTINUE
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL UPCUSR(WS,1)
      ENDIF
C
C
C  TEST FOR CORRECT CELL NUMBER AT COLLISION POINT
C  KILL PARTICLE, IF TOO LARGE ROUND OFF ERRORS DURING
C  PARTICLE TRACING
C
      IF (NLTEST) CALL CLLTST(*997)
      IF (ZT.GT.0.D0) THEN
C
C  FLIGHT WITH PARALLEL VELOCITY VELPAR (CM/SEC)
C  PARALLEL DISTANCE ZT (CM)
C  ENERGY RELAXATION CONSTANT TAUE
C
        DUR=ZT/VELPAR
        ENEW=E0*EXP(-DUR/TAUE)+1.5*TIIN(1,NCELL)*(1.-EXP(-DUR/TAUE))
        VNEW=RSQDVI(IOLD)*SQRT(ENEW)
C
C  UPDATE ESTIMATORS EIIO,EIPL
        EIIO(NCELLT)=EIIO(NCELLT)+WEIGHT*(ENEW-E0)
        EIPL(NCELLT)=EIPL(NCELLT)-WEIGHT*(ENEW-E0)
C
        E0=ENEW
        VEL=VNEW
        RETURN 2
      ENDIF
C
C  POST COLLISION ESTIMATOR
C
C     IF (NCLVI.GT.0) THEN
C       WS=WEIGHT/SIGTOT
C       CALL UPCUSR(WS,2)
C     ENDIF
      RETURN 2
C
997   CALL MASAGE ('ERROR IN FPKCOL,  DETECTED IN SUBR. CLLTST    ')
      CALL MASAGE ('PARTICLE IS KILLED                            ')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST
      IF (NLTRC) CALL CHCTRC(X0,Y0,Z0,16,18)
      GOTO 999
C
999   PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      WEIGHT=0.
      CALL LEER(1)
      RETURN 2
      END
