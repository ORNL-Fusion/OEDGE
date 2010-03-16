C  3D HISTOGRAM PLOT
C
      SUBROUTINE PL3DPG (ARR,IBLD,ICURV,
     .                   IX,IY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,LOGL,
     .                   ZMA,ZMI,W1,W2,
     .                   HEAD,RUNID,TXHEAD,TRC)

      USE PRECISION
      USE PARMMOD
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
C
      INTEGER, PARAMETER :: LAR=46*128*128

      REAL(DP), INTENT(IN) :: XX(*), YY(*)
      REAL(DP), INTENT(INOUT) :: ARR(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI, W1, W2
      INTEGER, INTENT(IN) :: IBLD, ICURV, IX, IY
      LOGICAL, INTENT(IN) :: LOGL,TRC
      CHARACTER(72), INTENT(IN) :: TEXT1, HEAD, RUNID, TXHEAD
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3

      REAL(DP) :: REMIN, REMAX, XMT, DX, YMI, YMA, RMI, RMA, AAR, XMINN,
     .          XMAXN, YMINN, YMAXN, YMT, XMI, XMA
!pb      REAL(SP) :: AR(LAR), EXT(3,3), VALU(3,2)
!pb      REAL(SP) :: XYZ(3,128,128)
!pb      real(sp) :: yh
      REAL :: AR(LAR), EXT(3,3), VALU(3,2)
      REAL :: XYZ(3,128,128)
      real :: yh
      INTEGER :: IR, IPX, IPY, IER, IPAN, IPEN, K, I, J
      CHARACTER(17) :: CH
      CHARACTER(20) :: CHAXS(3)
C
      DO 1 I=1,3
      DO 1 J=1,128
      DO 1 K=1,128
1       XYZ(I,J,K)=-75.75E20
C
      IF (LEVGEO.LE.1.OR.LEVGEO.GT.3.OR..NOT.LPTOR3(IBLD)) THEN
        WRITE (iunout,*) 'PLOTOPTION NOT READY. RETURN FROM PL3DPG '
        RETURN
      ENDIF
C
      IPAN=1
      IPEN=NP2ND
C
C  SEARCH MINIMA AND MAXIMA OF INDEPENDENT VARIABLES X AND Y
C
      XMI=1.D60
      XMA=-1.D60
      YMI=1.D60
      YMA=-1.D60
      RMI=1.D60
      RMA=-1.D60
      DO 10 I=1,NR1ST
        DO 10 J=1,NPPLG
          DO 10 K=NPOINT(1,J),NPOINT(2,J)
            XMI=MIN(XMI,XPOL(I,K))
            XMA=MAX(XMA,XPOL(I,K))
            YMI=MIN(YMI,YPOL(I,K))
            YMA=MAX(YMA,YPOL(I,K))
C  SEARCH MINIMA AND MAXIMA OF DEPENDENT VARIABLE Z=ARR
            IF (K.LT.NPOINT(2,J).AND.I.LT.NR1ST) THEN
              IR=I+(K-1)*NR1ST
              IF (LOGL) THEN
                RMI=MIN(RMI,MAX(1.E-48_DP,ARR(IR)))
                RMA=MAX(RMA,MAX(1.E-48_DP,ARR(IR)))
                ARR(IR)=LOG10(MAX(1.E-48_DP,ARR(IR)))
              ELSE
                RMI=MIN(RMI,ARR(IR))
                RMA=MAX(RMA,ARR(IR))
              ENDIF
            ENDIF
10    CONTINUE
C
      REMIN=ZMI
      REMAX=ZMA
      IF (LOGL) THEN
        IF (ZMI.EQ.666.) THEN
          REMIN=LOG10(RMI)
        ELSE
          REMIN=LOG10(MAX(1.E-48_DP,ZMI))
        ENDIF
        IF (ZMA.EQ.666.) THEN
          REMAX=LOG10(RMA)
        ELSE
          REMAX=LOG10(MAX(1.E-48_DP,ZMA))
        ENDIF
      ELSE
        IF (ZMI.EQ.666.) REMIN=RMI
        IF (ZMA.EQ.666.) REMAX=RMA
      ENDIF
C
      IF (TRC) THEN
        WRITE (iunout,*) 
     .    'PL3DPG:  ,IPAN,IPEN,XMI,XMA,YMI,YMA,REMIN,REMAX'
        WRITE (iunout,*) 
     .    '        ',IPAN,IPEN,XMI,XMA,YMI,YMA,REMIN,REMAX
      ENDIF
C
      IF (ABS(REMAX-REMIN)/MAX(REMAX,1.E-30_DP).LT.0.01)
     .    REMAX=REMIN+0.01*REMIN*SIGN(1._DP,REMIN)
C
C
      DX=MAX(ABS(XMA-XMI),ABS(YMA-YMI))*0.5
      XMT=0.5*(XMA+XMI)
      YMT=0.5*(YMA+YMI)
      XMINN=XMT-DX
      XMAXN=XMT+DX
      YMINN=YMT-DX
      YMAXN=YMT+DX
C
C     SET PLOT PARAMETER
C
C  DRDMPA(15) FOR EACH VALID PART, SEE BELOW
      CALL GRNXTB (1,'PL3DPG.F')
      CALL GRSCLC(10.,3.,34.,27.)
      CALL GRSCLV(2.,2.,26.,26.)
      CALL GR3DIM(LAR,IER)
C
      DO 30 K=1,NPPLG
        IPX=1
        DO 25 I=1,NR1ST-1
          IPY=1
          DO 20 J=NPOINT(1,K),NPOINT(2,K)-1
            IF (IPX+2.GT.128) GOTO 999
            IF (IPY+2.GT.128) GOTO 999
            XYZ(1,IPX+1,IPY+1)=XPOL(I,J)
            XYZ(1,IPX+1,IPY+2)=XPOL(I,J+1)
            XYZ(1,IPX+2,IPY+1)=XPOL(I+1,J)
            XYZ(1,IPX+2,IPY+2)=XPOL(I+1,J+1)
C
            XYZ(2,IPX+1,IPY+1)=YPOL(I,J)
            XYZ(2,IPX+1,IPY+2)=YPOL(I,J+1)
            XYZ(2,IPX+2,IPY+1)=YPOL(I+1,J)
            XYZ(2,IPX+2,IPY+2)=YPOL(I+1,J+1)
C
            IR=I+(J-1)*NR1ST
            AAR=ARR(IR)
            AAR=MAX(MIN(REMAX,AAR),REMIN)
            XYZ(3,IPX+1,IPY+1)=AAR
            XYZ(3,IPX+1,IPY+2)=AAR
            XYZ(3,IPX+2,IPY+1)=AAR
            XYZ(3,IPX+2,IPY+2)=AAR
            IPY=IPY+2
20        CONTINUE
          IPX=IPX+2
25      CONTINUE
C
        IF (IPX+1.GT.128) GOTO 999
        IF (IPY+1.GT.128) GOTO 999
        DO 26 I=2,IPY
          XYZ(1,1,I)=XYZ(1,2,I)
          XYZ(2,1,I)=XYZ(2,2,I)
          XYZ(3,1,I)=REMIN
          XYZ(1,IPX+1,I)=XYZ(1,IPX,I)
          XYZ(2,IPX+1,I)=XYZ(2,IPX,I)
          XYZ(3,IPX+1,I)=REMIN
26      CONTINUE
        DO 27 I=2,IPX
          XYZ(1,I,1)=XYZ(1,I,2)
          XYZ(2,I,1)=XYZ(2,I,2)
          XYZ(3,I,1)=REMIN
          XYZ(1,I,IPY+1)=XYZ(1,I,IPY)
          XYZ(2,I,IPY+1)=XYZ(2,I,IPY)
          XYZ(3,I,IPY+1)=REMIN
27      CONTINUE
        XYZ(1,1,1)=XYZ(1,2,2)
        XYZ(2,1,1)=XYZ(2,2,2)
        XYZ(3,1,1)=REMIN
        XYZ(1,IPX+1,1)=XYZ(1,IPX,2)
        XYZ(2,IPX+1,1)=XYZ(2,IPX,2)
        XYZ(3,IPX+1,1)=REMIN
        XYZ(1,IPX+1,IPY+1)=XYZ(1,IPX,IPY)
        XYZ(2,IPX+1,IPY+1)=XYZ(2,IPX,IPY)
        XYZ(3,IPX+1,IPY+1)=REMIN
        XYZ(1,1,IPY+1)=XYZ(1,2,IPY)
        XYZ(2,1,IPY+1)=XYZ(2,2,IPY)
        XYZ(3,1,IPY+1)=REMIN
        DO I=1,IPX+1
          DO J=1,IPY+1
            XYZ(1,I,J)=(XYZ(1,I,J)-XMINN)/(XMAXN-XMINN)
            XYZ(2,I,J)=(XYZ(2,I,J)-YMINN)/(YMAXN-YMINN)
            XYZ(3,I,J)=(XYZ(3,I,J)-REMIN)/(REMAX-REMIN)
          ENDDO
        ENDDO
        CALL GR3NET(AR,IER,128,XYZ,IPX+1,1,IPY+1,1,1,2)
30    CONTINUE
      CALL GR3EXT(AR,IER,EXT)
      VALU(1,1)=XMI
      VALU(1,2)=XMA
      VALU(2,1)=YMI
      VALU(2,2)=YMA
      VALU(3,1)=REMIN
      VALU(3,2)=REMAX
      CHAXS(1) = ' '
      CHAXS(2) = ' '
      CHAXS(3) = ' '
      CALL GR3AXS(AR,IER,EXT,VALU,CHAXS,.FALSE.,4,1)
      CALL GR3ROT(AR,IER,'Z',REAL(W1,KIND(1.E0)),
     .            'X',REAL(W2,KIND(1.E0)),'Y',0.0)
      CALL GR3PLO(AR,IER,'HID')
C
C     WRITE TEXT AND MEAN VALUE ONTO THE PLOT
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
      WRITE (CH,'(1P,E10.3)') RMA
      CALL GRTXT (1.,REAL(YH-2.5,KIND(1.E0)),10,CH)
      CALL GRTXT (1.,REAL(YH-3.,KIND(1.E0)),10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') RMI
      CALL GRTXT (1.,REAL(YH-3.5,KIND(1.E0)),10,CH)
C
      RETURN
999   CONTINUE
      WRITE (iunout,*) 'NOT ENOUGH STORAGE FOR 3D HISTOGRAM PLOT'
      WRITE (iunout,*) 'REDUCE PLOT AREA '
      WRITE (iunout,*) 'PLOT ABANDONNED'
      RETURN
      END
