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
     .           NZSAVE, INCZ, J1, IRSAVE, MTSAVE, LEARCA

      ZRADS=ZRAD
C
      IF (NLTRC) THEN
        CALL LEER(1)
        WRITE (6,*) 'TIMET: ZRAD,NRCELL,NTCELL,MTSURF '
        WRITE (6,*) '      ',ZRAD,NRCELL,NTCELL,MTSURF
        WRITE (6,*) 'INITIAL: X0,X01,Z01,PHI ',X0,X01,Z01,PHI/DEGRAD
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
            IF (NLTRC) WRITE (6,*) 'MTTEST,F,DZ ',MTTEST,F,DZ
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
                WRITE (6,*) 'ERROR FROM TIMET: '
                WRITE (6,*) 'RE-ENTRY THROUGH TOROIDAL SURFACE NOT '
                WRITE (6,*) 'READY FOR LEVGEO.GT.1 '
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
            WRITE (6,*) 'REENTRY FOUND, MTSURF,ZRAD = ',MTSURF,ZRAD
            WRITE (6,*) 'IRCELL ',IRCELL
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
            IF (NLTRC) WRITE (6,*) 'CHECK VALID RANGE: Z0TEST ',Z0TEST
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
            WRITE (6,*) 'NO REENTRY INTO TOROIDAL GRID FOUND '
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
          IF (NLTRC) WRITE (6,*) 'STOP AT MTSURF ',MTSURF
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
        WRITE (6,*) 'NPANU ',NPANU
        WRITE (6,*) 'X01,Z01 OUT OF RANGE IN TIMET'
        WRITE (6,*) X01,Z01,TTT
        WRITE (6,*) 'TRY TO KILL PARTICLE ASAP '
        ZRAD=-1._DP
        RETURN
      ENDIF
C
      Z001=Z01+ZRAD*VELZ
      X001=X01+ZRAD*VELX
      IF (ZRAD.LT.1.D30.AND.X01*X001.GT.0.D0) THEN
        ITT=IDINT(REAL(Z001/(X001*TANAL),KIND(1.D0)))
        IF (NLTRC) WRITE (6,*) 'TIMET 1 ',X01,Z01,X001,Z001,ITT
      ELSE
        TO=(Z01-X01*TANAL)/(TANAL*VELX-VELZ)
        XTO=X01+TO*VELX
        TU=(Z01+X01*TANAL)/(-TANAL*VELX-VELZ)
        XTU=X01+TU*VELX
        IF (NLTRC) WRITE (6,*) 'TU,TO ',TU,TO,XTU,XTO
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
        IF (NLTRC) WRITE (6,*) 'TIMET 2 ',X01,Z01,ITT
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
      IF (NLTRC) WRITE (6,*) 'FINAL 1: X01,Z01,PHI ',X01,Z01,PHI/DEGRAD
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
          IF (NLTRC) WRITE (6,*) 'TRY AGAIN IN TIMET'
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
        WRITE (6,*) 'FINAL 2: X01,Z01,PHI ',X01,Z01,PHI/DEGRAD
        WRITE (6,*) 'ZRAD,ISTS,ISRFCL,IPERID ',ZRAD,ISTS,ISRFCL,IPERID
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
          IF (NLTRC) WRITE (6,*) 'TRY AGAIN IN TIMET'
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
        WRITE (6,*) 'FINAL 3: X01,Z01,PHI ',X01,Z01,PHI/DEGRAD
        WRITE (6,*) 'ZRAD,ISTS,ISRFCL,IPERID ',ZRAD,ISTS,ISRFCL,IPERID
      ENDIF
      GOTO 5000
C
C  DISCRETE TOROIDAL APPROXIMATION FINISHED
C
C
5000  CONTINUE
      IF (NLTRC) WRITE (6,*) 'NCOUT= ',NCOUT
      DO 5100 J=1,NCOUT
        CLPD(J)=BLPD(J)
        NUPC(J)=(KUPC(J)-1)*NP2T3
        NCOUNT(J)=KUPC(J)
        IF (CLPD(J).LE.0..OR.KUPC(J).LE.0.OR.
     .      (KUPC(J).GE.NT3RD.AND.NLTOR)) THEN
          WRITE (6,*) 'ERROR DETECTED IN TIMET '
          WRITE (6,*) 'NPANU,J,BLPD,KUPC ',NPANU,J,BLPD(J),KUPC(J)
        ENDIF
        IF (NLTRC) THEN
          WRITE (6,*) 'TIMET: J,BLPD,NUPC,NCOUNT ',
     .                        J,BLPD(J),NUPC(J),NCOUNT(J)
        ENDIF
5100  CONTINUE
      IF (NLTRC) THEN
        WRITE (6,*) 'MTSURF,NLSRFZ,NINCZ,IRCELL,IPCELL '
        WRITE (6,*)  MTSURF,NLSRFZ,NINCZ,IRCELL,IPCELL
        IF (NLTOR) WRITE (6,*) 'INMP3I ',INMP3I(IRCELL,IPCELL,MTSURF)
      ENDIF
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
        WRITE (6,*) 'TRY TO KILL PARTICLE ASAP '
        SUM=-1.0D0
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
      CALL EXIT_OWN(1)
      END