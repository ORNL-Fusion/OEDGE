C
C
      SUBROUTINE RPSVEC (AORIG,BORIG,IBLD,ICURV,
     .                   IXX,IYY,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE PRODUCES A DATASET USED BY RAPS TO PRODUCE 
C  A VECTORFIELD PLOT
C
C  IT CALLS SUBR. CELINT, WHERE INTERPOLATION ONTO VERTICES IS PERFORMED
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CLOGAU
      USE CPLOT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CTRIG

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: AORIG(*), BORIG(*)
      REAL(DP), INTENT(IN) :: XX(*),YY(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI
      INTEGER, INTENT(IN) :: IBLD, ICURV, IXX, IYY
      LOGICAL, INTENT(IN) :: LOGL, TRC
      CHARACTER(72), INTENT(IN) :: TEXT1, HEAD, RUNID, TXHEAD
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3

      REAL(DP) :: XL, XM, BETRAG
      REAL(DP), ALLOCATABLE :: YWERT(:,:), YWERT1(:,:),
     .                       ZWERT(:,:), ZWERT1(:,:)
      INTEGER :: I, IP, IPART, IA, IB, IC, J, K, IT, LENCH, IERR, IR,
     .           NRAPS2, NVPLOT
      INTEGER :: ZUORD(NKNOT,0:20,2)
      REAL(SP) :: XY(800)
      REAL(SP) :: YH
      CHARACTER(17) :: CH
c
      data nvplot/0/
C
      WRITE (60,*) RUNID
      WRITE (60,*) TXHEAD
      WRITE (60,*) HEAD
      WRITE (60,*) TEXT1
      WRITE (60,*) TEXT2
      WRITE (60,*) TEXT3
      WRITE (60,*)
c
      nraps2=80
      nvplot=nvplot+1
      ch(1:7)='vecplot'
      if (nvplot.lt.10) then
        write (ch(8:8),'(i1)') nvplot
        lench=8
      else
        write (ch(8:9),'(i2)') nvplot
        lench=9
      endif
C
      NRAPS2=NRAPS2+nvplot
      nraps=nraps+1
      IRAPS=IRAPS+1
C
C  WRITE VALUE OF THE VECTOR TO RAPS-FILE IN ORDER TO HAVE A
C  SHADED PLOT
      OPEN (UNIT=NRAPS,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND NRAPS
C
      IF (LEVGEO.EQ.4) THEN
        ALLOCATE (YWERT1(NRAD,1))
        ALLOCATE (ZWERT1(NRAD,1))
        CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
        CALL CELINT(BORIG,ZWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
      ELSEIF (LEVGEO.LE.3) THEN
        ALLOCATE (YWERT(N1ST,N2ND+N3RD))
        ALLOCATE (ZWERT(N1ST,N2ND+N3RD))
        CALL CELINT(AORIG,YWERT,LOGL,IBLD,ICURV,N1ST,IERR)
        CALL CELINT(BORIG,ZWERT,LOGL,IBLD,ICURV,N1ST,IERR)
      ENDIF
      IF (IERR.GT.0) THEN
        IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
        IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
        IF (ALLOCATED(ZWERT)) DEALLOCATE (ZWERT)
        IF (ALLOCATED(ZWERT1)) DEALLOCATE (ZWERT1)
        RETURN
      END IF
C
      IF (LEVGEO.LE.2.AND.NLTOR) THEN
        DO 1100 IR=1,NR1ST
          DO 3100 IT=1,NT3RD
            IF (ZMI.NE.666.) YWERT(IR,IT)=MAX(YWERT(IR,IT),ZMI)
            IF (ZMA.NE.666.) YWERT(IR,IT)=MIN(YWERT(IR,IT),ZMA)
            IF (ZMI.NE.666.) ZWERT(IR,IT)=MAX(ZWERT(IR,IT),ZMI)
            IF (ZMA.NE.666.) ZWERT(IR,IT)=MIN(ZWERT(IR,IT),ZMA)
            BETRAG = SQRT(YWERT(IR,IT)**2+ZWERT(IR,IT)**2)
            WRITE (NRAPS,*) BETRAG
            YWERT(IR,IT)=YWERT(IR,IT)/(BETRAG+1.D-20)
            ZWERT(IR,IT)=ZWERT(IR,IT)/(BETRAG+1.D-20)
3100      CONTINUE
1100    CONTINUE
C
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.NLPOL) THEN
        DO 10 IR=1,NR1ST
          DO 20 IPART=1,NPPLG
            DO 30 IP=NPOINT(1,IPART),NPOINT(2,IPART)
              IF (ZMI.NE.666.) YWERT(IR,IP)=MAX(YWERT(IR,IP),ZMI)
              IF (ZMA.NE.666.) YWERT(IR,IP)=MIN(YWERT(IR,IP),ZMA)
              IF (ZMI.NE.666.) ZWERT(IR,IP)=MAX(ZWERT(IR,IP),ZMI)
              IF (ZMA.NE.666.) ZWERT(IR,IP)=MIN(ZWERT(IR,IP),ZMA)
              BETRAG = SQRT(YWERT(IR,IP)**2+ZWERT(IR,IP)**2)
              WRITE (NRAPS,*) BETRAG
              YWERT(IR,IP)=YWERT(IR,IP)/(BETRAG+1.D-20)
              ZWERT(IR,IP)=ZWERT(IR,IP)/(BETRAG+1.D-20)
30          CONTINUE
20        CONTINUE
10      CONTINUE
C
      ELSEIF (LEVGEO.EQ.4) THEN
        DO 60 I=1,NRKNOT
          IF (ZMI.NE.666.) YWERT1(I,1)=MAX(YWERT1(I,1),ZMI)
          IF (ZMA.NE.666.) YWERT1(I,1)=MIN(YWERT1(I,1),ZMA)
          IF (ZMI.NE.666.) ZWERT1(I,1)=MAX(ZWERT1(I,1),ZMI)
          IF (ZMA.NE.666.) ZWERT1(I,1)=MIN(ZWERT1(I,1),ZMA)
          BETRAG = SQRT(YWERT1(I,1)**2+ZWERT1(I,1)**2)
          WRITE (NRAPS,*) BETRAG
          YWERT1(I,1)=YWERT1(I,1)/(BETRAG+1.D-20)
          ZWERT1(I,1)=ZWERT1(I,1)/(BETRAG+1.D-20)
60      CONTINUE
C
      ELSE
        WRITE (6,*) 'UNWRITTEN OPTION IN RPSVEC: PLOT ABANDONNED '
        IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
        IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
        IF (ALLOCATED(ZWERT)) DEALLOCATE (ZWERT)
        IF (ALLOCATED(ZWERT1)) DEALLOCATE (ZWERT1)
        RETURN
      ENDIF
C
      CLOSE (UNIT=NRAPS)
C
C  WRITE VECTOR-COMPONENTS TO RAPS-FILE IN ORDER TO FORM A VECTORPLOT
C
C
      OPEN (UNIT=NRAPS2,file=ch(1:lench),
     .                  ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND NRAPS2

      WRITE(NRAPS2,'(1X,A5,8X,A4,50(11X,I1))') '-1111',
     .'PFEI',1,3,1,1,1
C
      IF (LEVGEO.LE.2.AND.NLTOR) THEN
        WRITE (6,*) 'UNWRITTEN OPTION IN RPSVEC: PLOT ABANDONNED '

      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.NLPOL) THEN
C
C       IR+1,IP+1          IR,IP+1        IR-1,IP+1_
C           +----------------+----------------+
C           |        5       |        4       |
C           | 6              |               3|
C           |                |                |
C       IR+1,IP ---------- IR,IP -------- IR-1,IP
C           |                |                |
C           | 7              |              2 |
C           |        8       |        1       |
C           + ---------------+----------------+
C       IR+1,IP-1          IR,IP-1        IR-1,IP-1
C
        DO IR=1,NR1ST
          DO IPART=1,NPPLG
            IPLOOP: DO IP=NPOINT(1,IPART),NPOINT(2,IPART)
C  RIGHT HEMISPHERE
              IF (IR.GT.1) THEN
C  LOWER RIGHT QUADRANT
                IF (IP.GT.NPOINT(1,IPART)) THEN
C  SIDE 1
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR,IP-1),YPOL(IR,IP-1),
     .                         XPOL(IR-1,IP-1),YPOL(IR-1,IP-1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
C  SIDE 2
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR-1,IP-1),YPOL(IR-1,IP-1),
     .                         XPOL(IR-1,IP),YPOL(IR-1,IP),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
                ENDIF
C  UPPER RIGHT QUADRANT
                IF (IP.LT.NPOINT(2,IPART)) THEN
C  SIDE 3
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR-1,IP),YPOL(IR-1,IP),
     .                         XPOL(IR-1,IP+1),YPOL(IR-1,IP+1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
C  SIDE 4
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR,IP+1),YPOL(IR,IP+1),
     .                         XPOL(IR-1,IP+1),YPOL(IR-1,IP+1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
                ENDIF
              ENDIF
C  LEFT HEMISPHERE
              IF (IR.LT.NR1ST) THEN
C  UPPER LEFT QUADRANT
                IF (IP.LT.NPOINT(2,IPART)) THEN
C  SIDE 5
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR+1,IP+1),YPOL(IR+1,IP+1),
     .                         XPOL(IR,IP+1),YPOL(IR,IP+1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
C  SIDE 6
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR+1,IP),YPOL(IR+1,IP),
     .                         XPOL(IR+1,IP+1),YPOL(IR+1,IP+1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
                ENDIF
C  LOWER LEFT QUADRANT
                IF (IP.GT.NPOINT(1,IPART)) THEN
C  SIDE 7
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR+1,IP),YPOL(IR+1,IP),
     .                         XPOL(IR+1,IP-1),YPOL(IR+1,IP-1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
C  SIDE 8
                  CALL MUELAM (XPOL(IR,IP),YPOL(IR,IP),
     .                         YWERT(IR,IP),ZWERT(IR,IP),
     .                         XPOL(IR+1,IP-1),YPOL(IR+1,IP-1),
     .                         XPOL(IR,IP-1),YPOL(IR,IP-1),XL,XM)
                  if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) GOTO 100
                ENDIF
              ENDIF
              CYCLE IPLOOP
100           CONTINUE
              ywert(IR,IP)=ywert(IR,IP)*xl*0.9
              zwert(IR,IP)=zwert(IR,IP)*xl*0.9
            ENDDO IPLOOP
          ENDDO
        ENDDO
        I=0
        DO IR=1,NR1ST
          DO IPART=1,NPPLG
            DO IP=NPOINT(1,IPART),NPOINT(2,IPART)
              I=I+1
              BETRAG=SQRT(YWERT(IR,IP)**2+ZWERT(IR,IP)**2)
              IF (BETRAG .GT. 1.E-5)
     .        WRITE(nraps2,'(I6,1P,5E12.4)')
     .             I,YWERT(IR,IP),zwert(IR,IP),0.,0.,0.
            enddo
          enddo
        enddo

      ELSEIF (LEVGEO.EQ.4) THEN
        DO 41 I=1,NRKNOT
          DO 51 J=0,20
            DO 51 K=1,2
            ZUORD(I,J,K) = 0
51        CONTINUE
41      CONTINUE
        DO 40 J=1,NTRII
          DO 50 I=1,3
            ZUORD(NECKE(I,J),0,1) = ZUORD(NECKE(I,J),0,1) + 1
c zuord darf maximal 20 werden
            if (zuord(necke(i,j),0,1).gt.20) then
              write (6,*) 'error in rpsvec: zuord'
              call exit_own(1)
            endif
            ZUORD(NECKE(I,J),ZUORD(NECKE(I,J),0,1),1) = J
            ZUORD(NECKE(I,J),ZUORD(NECKE(I,J),0,1),2) = I
50        CONTINUE
40      CONTINUE
        DO 61 I=1,NRKNOT
          IF (ZUORD(I,0,1).LT.1) THEN
            WRITE (6,*) 'ERROR IN RPSVEC: POINT ',I,' NOT IN MESH'
          ENDIF
          DO 70 J=1,ZUORD(I,0,1)
*  K: NUMBER OF TRIANGLE CONTAINING NODE I
            K=ZUORD(I,J,1)
*  IA: NUMBER OF NODE I IN TRIANGLE K
            ia=ZUORD(I,J,2)
c  ib, ic are the two other nodes in triangle k
            ib=ia+1
            if (ib.eq.4) ib=1
            ic=ib+1
            if (ic.eq.4) ic=1
c
            call muelam (xtrian(i),ytrian(i),ywert1(i,1),zwert1(i,1),
     .                   xtrian(necke(ib,k)),ytrian(necke(ib,k)),
     .                   xtrian(necke(ic,k)),ytrian(necke(ic,k)),xl,xm)
            if (xm >= 0. .and. xm <= 1. .AND. XL > 0.) then
              ywert1(i,1)=ywert1(i,1)*xl*0.9
              zwert1(i,1)=zwert1(i,1)*xl*0.9
              goto 61
            endif
70        CONTINUE
c   no intersection found.
c   point I is on a boundary, and the vector is pointing outside
c   the computational volume. don't plot it.
          YWERT1(I,1) = 0.
          YWERT1(I,1) = 0.
61      continue

        DO I=1,NRKNOT
          BETRAG=SQRT(YWERT1(I,1)**2+ZWERT1(I,1)**2)
          IF (BETRAG .GT. 1.D-5)
     .    WRITE(nraps2,'(I6,1P,5E12.4)') I,YWERT1(I,1),ZWERT1(I,1),
     .                                   0.,0.,0.
        enddo
      ENDIF
C
      WRITE(NRAPS2,'(1X,A5,8X,A3,50(11X,I1))') '-9999',
     .           'FIN',0,0,0
      CLOSE (UNIT=NRAPS2)
C
      IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
      IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
      IF (ALLOCATED(ZWERT)) DEALLOCATE (ZWERT)
      IF (ALLOCATED(ZWERT1)) DEALLOCATE (ZWERT1)
      RETURN
      END
