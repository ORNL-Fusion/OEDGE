C
      SUBROUTINE RPSCOL (AORIG,IBLD,ICURV,
     .                   I1,I2,XX,YY,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE PRODUCES A PLOTFILE FOR THE RAPS GRAPHICS SYSTEM
C
C  IT CALLS SUBR. CELINT, WHERE INTERPOLATION ONTO VERTICES IS PERFORMED
C
C
C    ******M******
C    *     *     *
C    *  D  I  E  *
C    *     *     *   |
C    L**H**A**F**J   |
C    *     *     *
C    *  C  G  B  *   IP
C    *     *     *
C    ******K******
C
C        <--- IR
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
      USE CCONA

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: AORIG(*)
      REAL(DP), INTENT(IN) :: XX(*), YY(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI
      INTEGER, INTENT(IN) :: IBLD, ICURV, I1, I2
      LOGICAL, INTENT(IN) :: LOGL, TRC
      CHARACTER(72), INTENT(IN) :: TEXT1, HEAD, RUNID, TXHEAD
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3
      
      REAL(DP), ALLOCATABLE :: YWERT(:,:),ywert1(:,:), BORIG(:)
      INTEGER :: IR, IERR, IT, I, IPART, IP, ICASE, IRD
C
      WRITE (60,*) RUNID
      WRITE (60,*) TXHEAD
      WRITE (60,*) HEAD
      WRITE (60,*) TEXT1
      WRITE (60,*) TEXT2
      WRITE (60,*) TEXT3
      WRITE (60,*)
C
      NRAPS=NRAPS+1
      IRAPS=IRAPS+1
C
      OPEN (UNIT=NRAPS,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND NRAPS
C
      ICASE=0
      IF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
         ALLOCATE (YWERT1(NRAD,1))
         CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
         ICASE=4
      ELSEIF (LEVGEO.LE.3.AND.
     .       (LPPOL3(IBLD).OR.LPTOR3(IBLD).OR.LPRAD3(IBLD))) THEN
         ALLOCATE (YWERT(N1ST,N2ND+N3RD))
         CALL CELINT(AORIG,YWERT,LOGL,IBLD,ICURV,N1ST,IERR)
         ICASE=3
      ELSEIF (LEVGEO.EQ.1.AND.NLTRZ.AND.
     .        NLRAD.AND.NLPOL.AND.NLTOR.AND..NOT.
     .       (LPPOL3(IBLD).OR.LPTOR3(IBLD).OR.LPRAD3(IBLD))) THEN
         ALLOCATE (YWERT1(NRAD,1))
         CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
         ICASE=1
      ELSEIF ((LEVGEO.EQ.5).AND..NOT.LRPSCUT) THEN
         ALLOCATE (YWERT1(NRAD,1))
         CALL CELINT(AORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
         ICASE=5
      ELSEIF ((LEVGEO.EQ.5).AND.LRPSCUT) THEN
         ALLOCATE (YWERT1(NRAD,1))
         ALLOCATE (BORIG(NRAD))
         CALL RPSCUT (AORIG,BORIG)
         CALL CELINT(BORIG,YWERT1,LOGL,IBLD,ICURV,NRAD,IERR)
         DEALLOCATE (BORIG)
         ICASE=4
      ENDIF

      IF (IERR.GT.0) THEN
        IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
        IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
        RETURN
      END IF
C
      IF (ICASE.EQ.1) THEN

         do ir=1,nr1st
            do ip=1,np2nd
               do it=1,nt3rd
                  IRD=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  IF (ZMI.NE.666.) YWERT1(IRD,1)=
     .                             MAX(YWERT1(IRD,1),ZMI)
                  IF (ZMA.NE.666.) YWERT1(IRD,1)=
     .                             MIN(YWERT1(IRD,1),ZMA)
                  IF (ABS(YWERT1(IRD,1)) < EPS30) YWERT1(IRD,1)=0._DP 
                  WRITE (NRAPS,*) YWERT1(IRD,1)
               enddo
            enddo
         enddo


      ELSEIF (LEVGEO.LE.2.AND.LPPOL3(IBLD)) THEN
        LPPOLR=.TRUE.
        IPPOLR=MAX(1,IPROJ3(IBLD,ICURV))
        IF (.NOT.NLRAD.OR..NOT.NLTOR) THEN
          WRITE (6,*) 'ERROR IN RPSCOL. LPPOL3? '
          IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
          IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
          RETURN
        ENDIF
        DO 1100 IR=1,NR1ST
          DO 3100 IT=1,NT3RD
            IF (ZMI.NE.666.) YWERT(IR,IT)=MAX(YWERT(IR,IT),ZMI)
            IF (ZMA.NE.666.) YWERT(IR,IT)=MIN(YWERT(IR,IT),ZMA)
            IF (ABS(YWERT(IR,IT)) < EPS30) YWERT(IR,IT)=0._DP 
            WRITE (NRAPS,*) YWERT(IR,IT)
3100      CONTINUE
1100    CONTINUE
C
      ELSEIF (LEVGEO.LE.2.AND.LPTOR3(IBLD)) THEN
        LPTORR=.TRUE.
        IPTORR=MAX(1,IPROJ3(IBLD,ICURV))
        IF (.NOT.NLRAD.OR..NOT.NLPOL) THEN
          WRITE (6,*) 'ERROR IN RPSCOL. LPTOR3? '
          IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
          IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
          RETURN
        ENDIF
        DO 1 IR=1,NR1ST
          DO 3 IP=1,NP2ND
            IF (ZMI.NE.666.) YWERT(IR,IP)=MAX(YWERT(IR,IP),ZMI)
            IF (ZMA.NE.666.) YWERT(IR,IP)=MIN(YWERT(IR,IP),ZMA)
            IF (ABS(YWERT(IR,IP)) < EPS30) YWERT(IR,IP)=0._DP 
            WRITE (NRAPS,*) YWERT(IR,IP)
3         CONTINUE
1       CONTINUE
C
      ELSEIF (LEVGEO.EQ.3.AND.LPTOR3(IBLD)) THEN
        LPTORR=.TRUE.
        IPTORR=MAX(1,IPROJ3(IBLD,ICURV))
        IF (.NOT.NLRAD.OR..NOT.NLPOL) THEN
          WRITE (6,*) 'ERROR IN RPSCOL. LPTOR3? '
          IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
          IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)
          RETURN
        ENDIF
        DO 10 IR=1,NR1ST
          DO 20 IPART=1,NPPLG
            DO 30 IP=NPOINT(1,IPART),NPOINT(2,IPART)
              IF (ZMI.NE.666.) YWERT(IR,IP)=MAX(YWERT(IR,IP),ZMI)
              IF (ZMA.NE.666.) YWERT(IR,IP)=MIN(YWERT(IR,IP),ZMA)
              IF (ABS(YWERT(IR,IP)) < EPS30) YWERT(IR,IP)=0._DP 
              WRITE (NRAPS,*) YWERT(IR,IP)
30          CONTINUE
20        CONTINUE
10      CONTINUE
C
C  icase=4
!PB      ELSEIF (LEVGEO.EQ.4.AND.LPTOR3(IBLD)) THEN
      ELSEIF (ICASE == 4) THEN
        LPTORR=.TRUE.
        IPTORR=MAX(1,IPROJ3(IBLD,ICURV))
        DO 60 I=1,NRKNOT
          IF (ZMI.NE.666.) YWERT1(I,1)=MAX(YWERT1(I,1),ZMI)
          IF (ZMA.NE.666.) YWERT1(I,1)=MIN(YWERT1(I,1),ZMA)
          IF (ABS(YWERT1(I,1)) < EPS30) YWERT1(I,1)=0._DP 
          WRITE(NRAPS,*) YWERT1(I,1)
60      CONTINUE
C
C  icase=5
      ELSEIF ((LEVGEO.EQ.5).AND..NOT.LRPSCUT) THEN
        DO I=1,NCOORD
          IF (ZMI.NE.666.) YWERT1(I,1)=MAX(YWERT1(I,1),ZMI)
          IF (ZMA.NE.666.) YWERT1(I,1)=MIN(YWERT1(I,1),ZMA)
          IF (ABS(YWERT1(I,1)) < EPS30) YWERT1(I,1)=0._DP 
          WRITE(NRAPS,*) YWERT1(I,1)
        enddo
C
      ELSE
        WRITE (6,*) 'UNWRITTEN OPTION IN RPSCOL: PLOT ABANDONNED '
      ENDIF
C
      CLOSE (UNIT=NRAPS)
C
      IF (ALLOCATED(YWERT)) DEALLOCATE (YWERT)
      IF (ALLOCATED(YWERT1)) DEALLOCATE (YWERT1)

      RETURN
      END



