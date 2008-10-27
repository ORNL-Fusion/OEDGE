

      SUBROUTINE ISOLNE (AORIG,IBLD,ICURV,
     .                   IXX,IYY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE CARRIES OUT A CONTOUR PLOT
C
C  IT CALLS SUBR. CELINT, WHERE INTERPOLATION ONTO VERTICES IS PERFORMED
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CLOGAU
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTRIG

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: AORIG(*)
      REAL(DP), INTENT(IN) :: XX(*), YY(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI
      INTEGER, INTENT(IN) :: IBLD, ICURV, IXX, IYY
      LOGICAL, INTENT(IN) :: LOGL, TRC
      CHARACTER(72), INTENT(IN) :: TEXT1, HEAD, RUNID, TXHEAD
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3

      REAL(DP) :: ZINT, SCLFCX, SCLFCY, RAMIN, RAMAX, FAK, CM, DX, DY,
     .          A1, A2, A3, A4, ACMIN, ACMAX, DA, ACONT, RMI, XMIN,
     .          XMAX, YMIN, YMAX, RMA, X1, X2, AA1, AA2
      REAL(DP), ALLOCATABLE :: A(:,:),AA(:,:)
      REAL(SP) :: XY(8000)
      REAL(SP) :: YH
      INTEGER :: NP, ITR, NP1, NP2, ICOLOR, IS, IC, IISO, NISO, IERR,
     .           IPART, IT, IR, IP, I
      CHARACTER(17) :: CH
C
      ZINT(X1,X2,AA1,AA2,A1)=X1+(A1-AA1)/(AA2-AA1+1.D-30)*(X2-X1)
C
C     PLOT 18 CONTOURS, WITH 6 DIFFERENT COLOURS
      NISO=18
      IISO=3
C
C  SEARCH FOR MAXIMUM AND MINIMUM RMI AND RMA
C
      RMI=1.D60
      RMA=-1.D60
C
      IF (LEVGEO .LE. 2.AND.LPTOR3(IBLD)) THEN
        IT=1
        IF (NLTOR) IT=IPROJ3(IBLD,ICURV)
        IF (IT.LE.0.OR.IT.GT.NT3RD) IT=1
        DO 21 IR=1,IXX-1
          DO 21 IP=1,IYY-1
            I=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            RMI=MIN(RMI,AORIG(I))
            RMA=MAX(RMA,AORIG(I))
21      CONTINUE
      ELSEIF (LEVGEO .LE. 2.AND.LPPOL3(IBLD)) THEN
        IP=1
        IF (NLPOL) IP=IPROJ3(IBLD,ICURV)
        IF (IP.LE.0.OR.IP.GT.NP2ND) IP=1
        DO 23 IR=1,IXX-1
          DO 23 IT=1,IYY-1
            I=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            RMI=MIN(RMI,AORIG(I))
            RMA=MAX(RMA,AORIG(I))
23      CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.LPTOR3(IBLD)) THEN
        IT=1
        IF (NLTOR) IT=IPROJ3(IBLD,ICURV)
        IF (IT.LE.0.OR.IT.GT.NT3RD) IT=1
        DO 20 IR=1,NR1ST-1
        DO 20 IPART=1,NPPLG
          DO 20 IP=NPOINT(1,IPART),NPOINT(2,IPART)-1
            I=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            RMI=MIN(RMI,AORIG(I))
            RMA=MAX(RMA,AORIG(I))
20      CONTINUE
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
        DO 22 I=1,NTRII
          RMI=MIN(RMI,AORIG(I))
          RMA=MAX(RMA,AORIG(I))
22      CONTINUE
      ELSE
        WRITE (6,*) 'MISSING OPTION IN ISOLNE: RMI,RMA '
        WRITE (6,*) 'PLOT ABANDONNED '
        RETURN
      ENDIF
C
C  INTERPOLATE: AORIG --> A
C
      IF (LEVGEO.EQ.4) THEN
        ALLOCATE (AA(NRAD,1))
        CALL CELINT(AORIG,AA,LOGL,IBLD,ICURV,NRAD,IERR)
      ELSEIF (LEVGEO.LE.3) THEN
         ALLOCATE (A(N1ST,MAX(N2ND,N3RD)))
        CALL CELINT(AORIG,A,LOGL,IBLD,ICURV,N1ST,IERR)
      ENDIF
      IF (IERR.GT.0) RETURN
C
C  SEARCH FOR XMIN,XMAX,YMIN,YMAX
C
      XMIN=1.D60
      YMIN=1.D60
      XMAX=-1.D60
      YMAX=-1.D60
C
      IF ((LEVGEO.LE.2).AND.LPPOL3(IBLD)) THEN
        XMIN = RHOSRF(1)
        XMAX = RHOSRF(NR1ST)
        YMIN = ZSURF(1)
        YMAX = ZSURF(NT3RD)
      ELSEIF (LEVGEO.EQ.1.AND.LPTOR3(IBLD)) THEN
        XMIN = RHOSRF(1)
        XMAX = RHOSRF(NR1ST)
        YMIN = PSURF(1)
        YMAX = PSURF(NP2ND)
      ELSEIF (LEVGEO.EQ.2.AND.LPTOR3(IBLD)) THEN
C  SUFFICIENT TO SEARCH ON OUTERMOST RADIAL SURFACE (BECAUSE: CONVEX)
        DO 5 IP = 1,NP2ND
          XMIN = MIN(XMIN,XPOL(NR1ST,IP))
          XMAX = MAX(XMAX,XPOL(NR1ST,IP))
          YMIN = MIN(YMIN,YPOL(NR1ST,IP))
          YMAX = MAX(YMAX,YPOL(NR1ST,IP))
5       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.LPTOR3(IBLD)) THEN
C  SEARCH ON WHOLE MESH
        DO 1 IR=1,NR1ST
          NP1=NPOINT(1,1)
          NP2=NPOINT(2,NPPLG)
          XMIN=MIN(XMIN,XPOL(IR,NP1),XPOL(IR,NP2))
          YMIN=MIN(YMIN,YPOL(IR,NP1),YPOL(IR,NP2))
          XMAX=MAX(XMAX,XPOL(IR,NP1),XPOL(IR,NP2))
          YMAX=MAX(YMAX,YPOL(IR,NP1),YPOL(IR,NP2))
1       CONTINUE
C
        DO 4 IPART=1,NPPLG
          DO 2 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            XMIN=MIN(XMIN,XPOL(1,IP))
            YMIN=MIN(YMIN,YPOL(1,IP))
            XMAX=MAX(XMAX,XPOL(1,IP))
            YMAX=MAX(YMAX,YPOL(1,IP))
2         CONTINUE
          DO 3 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            XMIN=MIN(XMIN,XPOL(NR1ST,IP))
            YMIN=MIN(YMIN,YPOL(NR1ST,IP))
            XMAX=MAX(XMAX,XPOL(NR1ST,IP))
            YMAX=MAX(YMAX,YPOL(NR1ST,IP))
3         CONTINUE
4       CONTINUE
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
C  SEARCH ON WHOLE MESH
        DO I=1,NRKNOT
          XMIN=MIN(XMIN,XTRIAN(I))
          YMIN=MIN(YMIN,YTRIAN(I))
          XMAX=MAX(XMAX,XTRIAN(I))
          YMAX=MAX(YMAX,YTRIAN(I))
        ENDDO
      ELSE
        WRITE (6,*) 'MISSING OPTION IN ISOLNE: XMIN,XMAX,YMIN,YMAX '
        WRITE (6,*) 'PLOT ABANDONNED '
        RETURN
      ENDIF
C
C
      CM=20.
      DX=(XMAX-XMIN)*FCABS1(IBLD)
      DY=(YMAX-YMIN)*FCABS2(IBLD)
      FAK=CM/MAX(DX,DY)
C
C  PLOT FRAME
C
      CALL GRNXTB (1)
      CALL GRSCLC (10.,4.,REAL(10.+DX*FAK,KIND(1.E0)),
     .                    REAL(4.+DY*FAK,KIND(1.E0)))
      CALL GRSCLV (REAL(XMIN,KIND(1.E0)),REAL(YMIN,KIND(1.E0)),
     .             REAL(XMAX,KIND(1.E0)),REAL(YMAX,KIND(1.E0)))
      CALL GRAXS (7,'X=3,Y=3',6,'R (CM)',6,'Z (CM)')
C
C  SCALE FACTORS: USER CO-ORDINATES TO CM:
C  X-DIRECTION:
      SCLFCX=((10.+DX*FAK)-10.)/(XMAX-XMIN)
C  Y-DIRECTION:
      SCLFCY=((4.+DY*FAK)-4.)/(YMAX-YMIN)
C
C  PLOT BOUNDARY OF MESH
C
      IF (LEVGEO.EQ.1) THEN
        CALL GRJMP(real(XMIN,KIND(1.E0)),real(YMIN,KIND(1.E0)))
        CALL GRDRW(real(XMIN,KIND(1.E0)),real(YMAX,KIND(1.E0)))
        CALL GRDRW(real(XMAX,KIND(1.E0)),real(YMAX,KIND(1.E0)))
        CALL GRDRW(real(XMAX,KIND(1.E0)),real(YMIN,KIND(1.E0)))
        CALL GRDRW(real(XMIN,KIND(1.E0)),real(YMIN,KIND(1.E0)))
      ELSEIF (LEVGEO.EQ.2.AND.LPPOL3(IBLD)) THEN
        CALL GRJMP(real(XMIN,KIND(1.E0)),real(YMIN,KIND(1.E0)))
        CALL GRDRW(real(XMIN,KIND(1.E0)),real(YMAX,KIND(1.E0)))
        CALL GRDRW(real(XMAX,KIND(1.E0)),real(YMAX,KIND(1.E0)))
        CALL GRDRW(real(XMAX,KIND(1.E0)),real(YMIN,KIND(1.E0)))
        CALL GRDRW(real(XMIN,KIND(1.E0)),real(YMIN,KIND(1.E0)))
      ELSEIF (LEVGEO.EQ.2.AND.LPTOR3(IBLD)) THEN
        DO 7 IR=1,NR1ST,NR1STM
          CALL GRJMP(real(XPOL(IR,1),KIND(1.E0)),
     .               real(YPOL(IR,1),KIND(1.E0)))
          DO 9 IP = 2,NP2ND
9           CALL GRDRW(real(XPOL(IR,IP),KIND(1.E0)),
     .                 real(YPOL(IR,IP),KIND(1.E0)))
7       CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.LPTOR3(IBLD)) THEN
        CALL GRJMP(real(XPOL(1,NPOINT(1,1)),KIND(1.E0)),
     .             real(YPOL(1,NPOINT(1,1)),KIND(1.E0)))
          DO 10 IR=2,NR1ST
10          CALL GRDRW (real(XPOL(IR,NPOINT(1,1)),KIND(1.E0)),
     .                  real(YPOL(IR,NPOINT(1,1)),KIND(1.E0)))
C
        CALL GRJMP (real(XPOL(1,NPOINT(2,NPPLG)),KIND(1.E0)),
     .              real(YPOL(1,NPOINT(2,NPPLG)),KIND(1.E0)))
        DO 11 IR=2,NR1ST
          NP=NPOINT(2,NPPLG)
11        CALL GRDRW (real(XPOL(IR,NP),KIND(1.E0)),
     .                real(YPOL(IR,NP),KIND(1.E0)))
        DO 15 I=1,NPPLG
          CALL GRJMP (real(XPOL(1,NPOINT(1,I)),KIND(1.E0)),
     .                real(YPOL(1,NPOINT(1,I)),KIND(1.E0)))
          DO 12 IP=NPOINT(1,I),NPOINT(2,I)
12          CALL GRDRW (real(XPOL(1,IP),KIND(1.E0)),
     .                  real(YPOL(1,IP),KIND(1.E0)))
          CALL GRJMP (real(XPOL(NR1ST,NPOINT(1,I)),KIND(1.E0)),
     .                real(YPOL(NR1ST,NPOINT(1,I)),KIND(1.E0)))
          DO 13 IP=NPOINT(1,I),NPOINT(2,I)
13          CALL GRDRW (real(XPOL(NR1ST,IP),KIND(1.E0)),
     .                  real(YPOL(NR1ST,IP),KIND(1.E0)))
15      CONTINUE
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
        DO ITR=1,NTRII
          IF (NCHBAR(1,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(1,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(1,ITR)),KIND(1.E0)))
            CALL GRDRW(REAL(XTRIAN(NECKE(2,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(2,ITR)),KIND(1.E0)))
          ENDIF
          IF (NCHBAR(2,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(3,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(3,ITR)),KIND(1.E0)))
            CALL GRDRW(REAL(XTRIAN(NECKE(2,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(2,ITR)),KIND(1.E0)))
          ENDIF
          IF (NCHBAR(3,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(1,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(1,ITR)),KIND(1.E0)))
            CALL GRDRW(REAL(XTRIAN(NECKE(3,ITR)),KIND(1.E0)),
     .                 REAL(YTRIAN(NECKE(3,ITR)),KIND(1.E0)))
          ENDIF
        ENDDO
      ELSE
        WRITE (6,*) 'MISSING OPTION IN ISOLNE: PLOT GRID BOUNDARY '
        WRITE (6,*) 'PLOT ABANDONNED '
        IF (ALLOCATED(A)) DEALLOCATE (A)
        IF (ALLOCATED(AA)) DEALLOCATE (AA)
        RETURN
      ENDIF
C
      IF (LOGL) THEN
        IF (ZMI.EQ.666.) THEN
          RAMIN=LOG10(MAX(1.E-48_DP,RMI))
        ELSE
          RAMIN=LOG10(MAX(1.E-48_DP,ZMI))
        ENDIF
        IF (ZMA.EQ.666.) THEN
          RAMAX=LOG10(MAX(1.E-48_DP,RMA))
        ELSE
          RAMAX=LOG10(MAX(1.E-48_DP,ZMA))
        ENDIF
      ELSE
        RAMIN=ZMI
        IF (ZMI.EQ.666.) RAMIN=RMI
        RAMAX=ZMA
        IF (ZMA.EQ.666.) RAMAX=RMA
      ENDIF
C
      IF (ABS(RAMAX-RAMIN)/MAX(RAMAX,1.E-30_DP).LT.0.01)
     .    RAMAX=RAMIN+0.01*RAMIN*SIGN(1._DP,RAMIN)
C
C
      DA=(RAMAX-RAMIN)/DBLE(NISO-1)
C
      ICOLOR=1
      IF ((LEVGEO.EQ.1).OR.(LEVGEO.EQ.2.AND.LPPOL3(IBLD))) THEN
        DO 1000 IS=1,NISO
          ACONT=RAMIN+(IS-1)*DA
          IC=0
          IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
          CALL GRNWPN(ICOLOR)
          DO 1100 IR=1,IXX-1
          DO 1100 IP=1,IYY-1
            A1=A(IR,IP)
            A2=A(IR+1,IP)
            A3=A(IR+1,IP+1)
            A4=A(IR,IP+1)
            ACMIN=MIN(A1,A2,A3,A4)
            ACMAX=MAX(A1,A2,A3,A4)
            IT=0
            IF (ACONT.GE.ACMIN.AND.ACONT.LE.ACMAX) THEN
              IF (ACONT.GE.MIN(A1,A2).AND.
     .            ACONT.LE.MAX(A1,A2)) THEN
                XY(IC+1)=ZINT(XX(IR),XX(IR+1),A1,A2,ACONT)
                XY(IC+2)=YY(IP)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A2,A3).AND.
     .            ACONT.LE.MAX(A2,A3)) THEN
                XY(IC+1)=XX(IR+1)
                XY(IC+2)=ZINT(YY(IP),YY(IP+1),A2,A3,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A3,A4).AND.
     .            ACONT.LE.MAX(A3,A4)) THEN
                XY(IC+1)=ZINT(XX(IR+1),XX(IR),A3,A4,ACONT)
                XY(IC+2)=YY(IP+1)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A4,A1).AND.
     .            ACONT.LE.MAX(A4,A1)) THEN
                XY(IC+1)=XX(IR)
                XY(IC+2)=ZINT(YY(IP+1),YY(IP),A4,A1,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (MOD(IT,2).EQ.1) THEN
C               WRITE (6,*) ' ACHTUNG !!!! '
C               WRITE (6,*) IT,' PUNKTE AUF DEM VIERECK GEFUNDEN '
C               WRITE (6,*) ' EIN ZUSAETZLICHER PUNKT EINGEGEBEN '
                XY(IC+1)=XY(IC-IT*2+1)
                XY(IC+2)=XY(IC-IT*2+2)
                IC=IC+2
              ENDIF
              IF (IC+8.GT.800) THEN
                CALL XYPLOT (XY,IC)
                IC=0
              ENDIF
            ENDIF
1100      CONTINUE
          IF (IC.GT.0) CALL XYPLOT (XY,IC)
1000    CONTINUE
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTOR3(IBLD)) THEN
        DO 100 IS=1,NISO
          ACONT=RAMIN+(IS-1)*DA
          IC=0
          IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
          CALL GRNWPN(ICOLOR)
          DO 110 IR=1,NR1ST-1
          DO 110 IPART=1,NPPLG
          DO 110 IP=NPOINT(1,IPART),NPOINT(2,IPART)-1
            A1=A(IR,IP)
            A2=A(IR+1,IP)
            A3=A(IR+1,IP+1)
            A4=A(IR,IP+1)
            ACMIN=MIN(A1,A2,A3,A4)
            ACMAX=MAX(A1,A2,A3,A4)
            IT=0
            IF (ACONT.GE.ACMIN.AND.ACONT.LE.ACMAX) THEN
              IF (ACONT.GE.MIN(A1,A2).AND.
     .            ACONT.LE.MAX(A1,A2)) THEN
                XY(IC+1)=ZINT(XPOL(IR,IP),XPOL(IR+1,IP),A1,A2,ACONT)
                XY(IC+2)=ZINT(YPOL(IR,IP),YPOL(IR+1,IP),A1,A2,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A2,A3).AND.
     .            ACONT.LE.MAX(A2,A3)) THEN
                XY(IC+1)=ZINT(XPOL(IR+1,IP),XPOL(IR+1,IP+1),A2,A3,ACONT)
                XY(IC+2)=ZINT(YPOL(IR+1,IP),YPOL(IR+1,IP+1),A2,A3,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A3,A4).AND.
     .            ACONT.LE.MAX(A3,A4)) THEN
                XY(IC+1)=ZINT(XPOL(IR+1,IP+1),XPOL(IR,IP+1),A3,A4,ACONT)
                XY(IC+2)=ZINT(YPOL(IR+1,IP+1),YPOL(IR,IP+1),A3,A4,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A4,A1).AND.
     .            ACONT.LE.MAX(A4,A1)) THEN
                XY(IC+1)=ZINT(XPOL(IR,IP+1),XPOL(IR,IP),A4,A1,ACONT)
                XY(IC+2)=ZINT(YPOL(IR,IP+1),YPOL(IR,IP),A4,A1,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (MOD(IT,2).EQ.1) THEN
C               WRITE (6,*) ' ACHTUNG !!!! '
C               WRITE (6,*) IT,' PUNKTE AUF DEM VIERECK GEFUNDEN '
C               WRITE (6,*) ' EIN ZUSAETZLICHER PUNKT EINGEGEBEN '
                XY(IC+1)=XY(IC-IT*2+1)
                XY(IC+2)=XY(IC-IT*2+2)
                IC=IC+2
              ENDIF
              IF (IC+8.GT.800) THEN
                CALL XYPLOT (XY,IC)
                IC=0
              ENDIF
            ENDIF
110       CONTINUE
          IF (IC.GT.0) CALL XYPLOT (XY,IC)
100     CONTINUE
      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
        DO IS=1,NISO
          ACONT=RAMIN+(IS-1)*DA
          IC=0
          IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
          CALL GRNWPN(ICOLOR)
          DO  ITR=1,NTRII
            A1=AA(NECKE(1,ITR),1)
            A2=AA(NECKE(2,ITR),1)
            A3=AA(NECKE(3,ITR),1)
            ACMIN=MIN(A1,A2,A3)
            ACMAX=MAX(A1,A2,A3)
            IT=0
            IF (ACONT.GE.ACMIN.AND.ACONT.LE.ACMAX) THEN
              IF (ACONT.GE.MIN(A1,A2).AND.
     .            ACONT.LE.MAX(A1,A2)) THEN
                XY(IC+1)=ZINT(XTRIAN(NECKE(1,ITR)),
     .                        XTRIAN(NECKE(2,ITR)),A1,A2,ACONT)
                XY(IC+2)=ZINT(YTRIAN(NECKE(1,ITR)),
     .                        YTRIAN(NECKE(2,ITR)),A1,A2,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A2,A3).AND.
     .            ACONT.LE.MAX(A2,A3)) THEN
                XY(IC+1)=ZINT(XTRIAN(NECKE(2,ITR)),
     .                        XTRIAN(NECKE(3,ITR)),A2,A3,ACONT)
                XY(IC+2)=ZINT(YTRIAN(NECKE(2,ITR)),
     .                        YTRIAN(NECKE(3,ITR)),A2,A3,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (ACONT.GE.MIN(A3,A1).AND.
     .            ACONT.LE.MAX(A3,A1)) THEN
                XY(IC+1)=ZINT(XTRIAN(NECKE(3,ITR)),
     .                        XTRIAN(NECKE(1,ITR)),A3,A1,ACONT)
                XY(IC+2)=ZINT(YTRIAN(NECKE(3,ITR)),
     .                        YTRIAN(NECKE(1,ITR)),A3,A1,ACONT)
                IT=IT+1
                IC=IC+2
              ENDIF
              IF (MOD(IT,2).EQ.1) THEN
C               WRITE (6,*) ' ACHTUNG !!!! '
C               WRITE (6,*) IT,' PUNKTE AUF DEM VIERECK GEFUNDEN '
C               WRITE (6,*) ' EIN ZUSAETZLICHER PUNKT EINGEGEBEN '
                XY(IC+1)=XY(IC-IT*2+1)
                XY(IC+2)=XY(IC-IT*2+2)
                IC=IC+2
              ENDIF
              IF (IC+8.GT.800) THEN
                CALL XYPLOT (XY,IC)
                IC=0
              ENDIF
            ENDIF
          ENDDO
          IF (IC.GT.0) CALL XYPLOT (XY,IC)
        ENDDO
      ELSE
        WRITE (6,*) 'MISSING OPTION IN ISOLNE: PLOT CONTOURS '
        WRITE (6,*) 'PLOT ABANDONNED '
        IF (ALLOCATED(A)) DEALLOCATE (A)
        IF (ALLOCATED(AA)) DEALLOCATE (AA)
        RETURN
      ENDIF
C
C
C     WRITE TEXT AND MEAN VALUE ONTO THE PLOT
C
      CALL GRNWPN (1)
      CALL GRSCLC (0.,0.,39.,28.)
      CALL GRSCLV (0.,0.,39.,28.)
      YH=27.5
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,RUNID)
      YH=26.75
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,HEAD)
      YH=26.00
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,REAL(YH,KIND(1.E0)),10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,REAL(YH-0.5,KIND(1.E0)),10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,REAL(YH-1.,KIND(1.E0)),10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,REAL(YH-2.,KIND(1.E0)),10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') RMA
      CALL GRTXT (1.,REAL(YH-2.5,KIND(1.E0)),10,CH)
      CALL GRTXT (1.,REAL(YH-3.,KIND(1.E0)),10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') RMI
      CALL GRTXT (1.,REAL(YH-3.5,KIND(1.E0)),10,CH)
C
      ICOLOR=1
      YH=YH-4.
      DO 200 IS=1,NISO
        ACONT=RAMIN+(IS-1)*DA
        IF (LOGL) ACONT=10.**ACONT
        IF (MOD(IS,IISO).EQ.1) ICOLOR=ICOLOR+1
        CALL GRNWPN(ICOLOR)
        YH=YH-0.5
        CALL GRJMP (1.,REAL(YH+0.25,KIND(1.E0)))
        CALL GRDRW (2.5,REAL(YH+0.25,KIND(1.E0)))
        CALL GRNWPN (1)
        WRITE (CH,'(1P,E10.3)') ACONT
        CALL GRTXT (3.,REAL(YH,KIND(1.E0)),10,CH)
200   CONTINUE
C
      IF (ALLOCATED(A)) DEALLOCATE (A)
      IF (ALLOCATED(AA)) DEALLOCATE (AA)

      RETURN
      END
