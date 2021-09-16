C
C
      SUBROUTINE PLOT3D (ARR,IBLD,ICURV,
     .                   NX,NY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,LOGL,
     .                   ZMA,ZMI,W1,W2,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  ON INPUT:  LOGL: USE LOG SCALE FOR ORDINATE
C             ARR(IJ),I=1,NX-1,J=1,NY-1 ,IJ=I+(J-1)*NR1ST
C                          ARRAY TO BE PLOTTED
C             XX(I),I=1,NX X-GRID BOUNDARIES
C             YY(J),J=1,NY Y GRID BOUNDARIES
C  ARR(IJ) IS THEN SET ONTO 2D ARRAY FALT(I,J)
C
C
C  LPOLAR: R-THETA CO-ORDINATES
C  LKARTH: X-Y     CO-ORDINATES
C
C  FALT-->XZY (3,...)
C
      USE PRECISION
      USE PARMMOD
      USE CPLOT
      USE CGRID
      USE CGEOM

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XX(*), YY(*)
      REAL(DP), INTENT(INOUT) :: ARR(*)
      REAL(DP), INTENT(IN) :: ZMI, ZMA, W1, W2
      INTEGER, INTENT(IN) :: NX, NY, IBLD, ICURV
      LOGICAL, INTENT(IN) :: LOGL,TRC
      CHARACTER(72), INTENT(IN) :: TEXT1
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3
      CHARACTER(72), INTENT(IN) :: HEAD, RUNID, TXHEAD

      REAL(DP) :: XMINN, XMAXN, DXXX, YMINN, YMAXN
!pb      REAL(SP) :: FALT(N1ST,N2NDPLG),
!pb     .          X(N1ST,N2NDPLG), Y(N1ST,N2NDPLG), Z2(N1ST,N2NDPLG)
!pb      REAL(SP) :: EXT(3,3), VALU(3,2), DCM, YH
!pb      REAL(SP) :: YHLF, XHLF, FMIN, FMAX, REMIN, REMAX,
!pb     .          XMIN, XMAX, YMIN, YMAX
!pb      REAL(SP), ALLOCATABLE :: AR(:)
      REAL :: FALT(N1ST,N2NDPLG),
     .        X(N1ST,N2NDPLG), Y(N1ST,N2NDPLG), Z2(N1ST,N2NDPLG)
      REAL :: EXT(3,3), VALU(3,2), DCM, YH
      REAL :: YHLF, XHLF, FMIN, FMAX, REMIN, REMAX,
     .        XMIN, XMAX, YMIN, YMAX
      REAL, ALLOCATABLE :: AR(:)
      INTEGER :: I, J, IXM, IYM, IZ, IZN, IER, LAR, IX, IY
      CHARACTER(17) :: CH
      CHARACTER(20) :: CHAXS(3)
C
C
      LAR=46*N1ST*N2NDPLG
      ALLOCATE (AR(LAR))
      ier=1
      IF (LEVGEO.LE.3.AND.LEVGEO.GT.1.AND.LPTOR3(IBLD)) THEN
C
        IXM=NX-1
        IYM=NY-1
        IZ=IYM*NR1ST
C
        DO 1 I=1,N1ST
          DO 1 J=1,N2NDPLG
            FALT(I,J)=-75.75E20
1       CONTINUE
C
        IF (LOGL) THEN
          DO 3 J=1,IZ
            ARR(J)=LOG10(MAX(1.E-48_DP,ARR(J)))
3         CONTINUE
        ENDIF
C
C  SET ONTO 2D ARRAY FOR PLOTTING
C
        DO 20 I=1,NX
          DO 20 J=1,IYM
            DXXX=ARR(I+(J-1)*NR1ST)
            FALT(I,J)=DXXX
            X(I,J)=XPOL(I,J)
            Y(I,J)=YPOL(I,J)
20      CONTINUE
C
        IXM=NX-1
        IYM=NY-1
C
C     SEARCH FOR MINIMUM AND MAXIMUM AND REPLACE, IF REQUIRED
C
        IZN=IXM*IYM
        FMIN=FALT(1,1)
        FMAX=FALT(1,1)
        XMIN=X(1,1)
        XMAX=X(1,1)
        YMIN=Y(1,1)
        YMAX=Y(1,1)
        DO 25 J=1,IXM
        DO 25 I=1,IYM
          FMIN=MIN(FMIN,FALT(J,I))
          FMAX=MAX(FMAX,FALT(J,I))
          XMIN=MIN(XMIN,X(J,I))
          XMAX=MAX(XMAX,X(J,I))
          YMIN=MIN(YMIN,Y(J,I))
          YMAX=MAX(YMAX,Y(J,I))
25      CONTINUE
C
        REMIN=ZMI
        REMAX=ZMA
        IF (LOGL) THEN
          REMIN=LOG10(MAX(1.E-48_DP,ZMI))
          REMAX=LOG10(MAX(1.E-48_DP,ZMA))
        ENDIF
        IF (ZMI.EQ.666.) REMIN=FMIN
        IF (ZMA.EQ.666.) REMAX=FMAX
        IF (REMIN.GE.REMAX) THEN
          REMIN=REMIN-1.
          REMAX=REMAX+1.
        ENDIF
C
        DO 30 J=1,IXM
        DO 30 I=1,IYM
          FALT(J,I)=MIN(REMAX,FALT(J,I))
          FALT(J,I)=MAX(REMIN,FALT(J,I))
30      CONTINUE
C
C  PLOT SMOOTH SURFACE
C
        DCM=MAX(ABS(XMAX-XMIN),ABS(YMAX-YMIN))
        XHLF=0.5*(XMAX+XMIN)
        YHLF=0.5*(YMAX+YMIN)
        XMINN=XHLF-0.5*DCM
        XMAXN=XHLF+0.5*DCM
        YMINN=YHLF-0.5*DCM
        YMAXN=YHLF+0.5*DCM
C
C     NORMIEREN DER WERTE
C
        DO IX=1,IXM
          DO IY=1,IYM
            X(IX,IY) = (X(IX,IY)-XMINN)/(XMAXN-XMINN)
            Y(IX,IY) = (Y(IX,IY)-YMINN)/(YMAXN-YMINN)
            FALT(IX,IY) = (FALT(IX,IY)-REMIN)/(REMAX-REMIN)
            Z2(IX,IY) = 0.
          ENDDO
        ENDDO
C
        CALL GRNXTB (1)
        CALL GRSCLC(10.,3.,34.,27.)
        CALL GRSCLV(2.,2.,26.,26.)
        CALL GR3DIM(LAR,IER)
        CALL GR3NT1(AR,IER,n1st,X,Y,Z2,IXM,1,IYM,1,1,1)
        CALL GR3NT1(AR,IER,n1st,X,Y,FALT,IXM,1,IYM,1,1,2)
        CALL GR3EXT(AR,IER,EXT)
        VALU(1,1)=XMIN
        VALU(1,2)=XMAX
        VALU(2,1)=YMIN
        VALU(2,2)=YMAX
        VALU(3,1)=REMIN
        VALU(3,2)=REMAX
        CHAXS(1) = ' '
        CHAXS(2) = ' '
        CHAXS(3) = ' '
        CALL GR3AXS(AR,IER,EXT,VALU,CHAXS,.FALSE.,4,1)
        CALL GR3ROT(AR,IER,'Z',REAL(W1,KIND(1.E0)),'X',
     .              REAL(W2,KIND(1.E0)),'Y',0.0)
        CALL GR3PLO(AR,IER,'HID')
C
      ELSE
        WRITE (6,*) 'INVALID OPTION IN PLOT3D  '
        WRITE (6,*) '3D-PLOT ABANDONED  '
        RETURN
      ENDIF
C
C     WRITE TEXT ONTO THE PLOT
C
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
      WRITE (CH,'(1P,E10.3)') FMAX
      CALL GRTXT (1.,REAL(YH-2.5,KIND(1.E0)),10,CH)
      CALL GRTXT (1.,REAL(YH-3.,KIND(1.E0)),10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') FMIN
      CALL GRTXT (1.,REAL(YH-3.5,KIND(1.E0)),10,CH)

      DEALLOCATE (AR)
C
      RETURN
      END
