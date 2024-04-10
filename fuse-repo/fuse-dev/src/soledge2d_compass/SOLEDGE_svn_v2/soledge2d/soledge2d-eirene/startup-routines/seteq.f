C
C
! 6.12.05  Bugfix: one index I replaced by J. Bug was relevant for switching of
!                  IGJUM2 in case of bit arrays
! 1.02.08  Bugfix: check of IGJUM1 replaced by check of IGJUM2
!                  Bug was relevant for switching of IGJUM2 in case of bit arrays
 
      SUBROUTINE EIRENE_SETEQ
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_CTRCEI
      USE EIRMOD_CLGIN
 
      IMPLICIT NONE
 
      REAL(DP) :: SG
      INTEGER :: NLBT, IEQ, J, IDIMP, K, I, IEQ1
      LOGICAL EIRENE_BITGET
C
C  FIRST: ADDITIONAL SURFACES
C
      DO 1 J=1,NLIMI
        IEQ=IABS(ILEQUI(J))
        IF (IEQ.EQ.0) GOTO 1
        IF (IGJUM0(J)+IGJUM0(IEQ) .NE. 0) GOTO 1
        SG=ISIGN(1,ILEQUI(J))
C   COEFFICIENTS OF SURFACE J ARE SET EQUAL TO THOSE OF SURFACE IEQ
        A0LM(J)=SG*A0LM(IEQ)
        A1LM(J)=SG*A1LM(IEQ)
        A2LM(J)=SG*A2LM(IEQ)
        A3LM(J)=SG*A3LM(IEQ)
        A4LM(J)=SG*A4LM(IEQ)
        A5LM(J)=SG*A5LM(IEQ)
        A6LM(J)=SG*A6LM(IEQ)
        A7LM(J)=SG*A7LM(IEQ)
        A8LM(J)=SG*A8LM(IEQ)
        A9LM(J)=SG*A9LM(IEQ)
        IF (JUMLIM(J).EQ.0) THEN
          ALM(J)=SG*ALM(IEQ)
          BLM(J)=SG*BLM(IEQ)
          CLM(J)=SG*CLM(IEQ)
        ELSE
          ALM(J)=ALM(IEQ)
          BLM(J)=BLM(IEQ)
          CLM(J)=CLM(IEQ)
        ENDIF
        IF (JUMLIM(J).NE.0) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM1(J,IEQ)=1
            IGJUM1(IEQ,J)=1
          ELSE
            CALL EIRENE_BITSET (IGJUM1,0,NLIMPS,J,IEQ,1,NBITS)
            CALL EIRENE_BITSET (IGJUM1,0,NLIMPS,IEQ,J,1,NBITS)
          END IF
        ELSE
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM2(J,IEQ)=1
            IGJUM2(IEQ,J)=1
          ELSE
            CALL EIRENE_BITSET (IGJUM2,0,NLIMPS,J,IEQ,1,NBITS)
            CALL EIRENE_BITSET (IGJUM2,0,NLIMPS,IEQ,J,1,NBITS)
          END IF
        ENDIF
C
        IF (TRCSUR) THEN
          WRITE (iunout,*) 'EQUATION OF SURFACE NO. ',J,
     .                     ' IS SET EQUAL TO '
          WRITE (iunout,*) 'EQUATION OF SURFACE NO. ',IEQ
          CALL EIRENE_LEER(1)
        ENDIF
1     CONTINUE
C
C  NEXT: NON DEFAULT STANDARD SURFACES
C  JUMLIM=0 FOR STANDARD SURFACES, ONLY LGJUM2 IS USED IN TIME-ROUTINES
C
      DO 10 J=NLIM+1,NLIM+NSTSI
        IF (ILEQUI(J).EQ.0.OR.IGJUM0(J).NE.0) GOTO 10
        IF (INUMP(J-NLIM,1).NE.0) IDIMP=1
        IF (INUMP(J-NLIM,2).NE.0) IDIMP=2
        IF (INUMP(J-NLIM,3).NE.0) IDIMP=3
C  NEGATIVE SIGN OF ILEQUI IS MEANINGLESS FOR STANDARD SURFACES
C  BECAUSE THEIR ORIENTATION CANNOT BE CHANGED
        IEQ1=IABS(ILEQUI(J))
C  FIND INDEX IEQ OF SURFACE IEQ1, TO BE EQUALLED WITH SURFACE J
C  SEARCH ONLY IN THE SAME (RADIAL, POLOIDAL OR TOROIDAL) GRID
        IEQ=0
        DO 11 I=1,NSTSI
          IF (INUMP(I,IDIMP).EQ.IEQ1) IEQ=NLIM+I
11      CONTINUE
        IF (NLIMPB >= NLIMPS) THEN
          IGJUM2(J,IEQ)=1
          IGJUM2(IEQ,J)=1
        ELSE
          CALL EIRENE_BITSET (IGJUM2,0,NLIMPS,J,IEQ,1,NBITS)
          CALL EIRENE_BITSET (IGJUM2,0,NLIMPS,IEQ,J,1,NBITS)
        END IF
        IF (TRCSUR) THEN
          WRITE (iunout,*) 'EQUATION OF SURFACE NO. ',J,
     .                     ' IS SET EQUAL TO '
          WRITE (iunout,*) 'EQUATION OF SURFACE NO. ',IEQ
          CALL EIRENE_LEER(1)
        ENDIF
10    CONTINUE
C
C  LGJUM1: COMMUTATIVE
C  LGJUM2: ASSOCIATIVE
      IF (NLIMPB >= NLIMPS) THEN
        DO J=1,NLIMI
          IF (SUM(IGJUM1(J,1:NLIMI)) > 1 ) THEN
            DO I=1,NLIMI
              IF (IGJUM1(J,I) .NE. 0) IGJUM1(I,J)=1
            ENDDO
          ENDIF
          IF (JUMLIM(J).EQ.0) THEN
            IF (SUM(IGJUM2(J,1:NLIMI)) > 1 ) THEN
              DO I=1,NLIMI
                DO K=1,NLIMI
                  IF (IGJUM2(I,J)*IGJUM2(J,K) .NE. 0) THEN
                    IGJUM2(I,K)=1
                    IGJUM2(K,I)=1
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ELSE
        NLBT = NLIMI/NBITS+1
        DO J=1,NLIMI
          IF (SUM(IGJUM1(J,1:NLBT)) > 1 ) THEN
            DO I=1,NLIMI
              IF (EIRENE_BITGET(IGJUM1,0,NLIMPS,J,I,NBITS))
     .          CALL EIRENE_BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
            END DO
          END IF
          IF (JUMLIM(J).EQ.0) THEN
            IF (SUM(IGJUM2(J,1:NLBT)) > 1 ) THEN
              DO I=1,NLIMI
                IF (EIRENE_BITGET(IGJUM2,0,NLIMPS,I,J,NBITS)) THEN
                  DO K=1,NLIMI
!pb                    IF (EIRENE_BITGET(IGJUM1,0,NLIMPS,J,K,NBITS)) THEN
                    IF (EIRENE_BITGET(IGJUM2,0,NLIMPS,J,K,NBITS)) THEN
                      CALL EIRENE_BITSET (IGJUM2,0,NLIMPS,I,K,1,NBITS)
                      CALL EIRENE_BITSET (IGJUM2,0,NLIMPS,K,I,1,NBITS)
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END IF
        END DO
      END IF
 
!pb apparently not needed as IGJUM2 refers to second order additional
!pb surfaces only
      IF (.FALSE.) THEN
      DO 200 J=NLIM+1,NLIM+NSTSI
        DO 200 I=NLIM+1,NLIM+NSTSI
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM1(J,I).NE.0) IGJUM1(I,J)=1
          DO K=NLIM+1,NLIM+NSTSI
            IF (IGJUM2(I,J)*IGJUM2(J,K) .NE. 0) THEN
              IGJUM2(I,K)=1
              IGJUM2(K,I)=1
            ENDIF
          END DO
        ELSE
          IF (EIRENE_BITGET(IGJUM1,0,NLIMPS,J,I,NBITS))
     .      CALL EIRENE_BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
          DO K=NLIM+1,NLIM+NSTSI
            IF (EIRENE_BITGET(IGJUM2,0,NLIMPS,I,J,NBITS) .AND.
     .          EIRENE_BITGET(IGJUM1,0,NLIMPS,J,K,NBITS)) THEN
              CALL EIRENE_BITSET (IGJUM2,0,NLIMPS,I,K,1,NBITS)
              CALL EIRENE_BITSET (IGJUM2,0,NLIMPS,K,I,1,NBITS)
            END IF
          END DO
        END IF
200   CONTINUE
      END IF
      RETURN
      END
