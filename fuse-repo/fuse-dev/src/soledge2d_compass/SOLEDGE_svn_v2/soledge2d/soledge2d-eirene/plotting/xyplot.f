C
C
      SUBROUTINE EIRENE_XYPLOT (XY,NC)
 
      USE EIRMOD_PRECISION
 
      IMPLICIT NONE
C
      REAL(SP), INTENT(INOUT) :: XY(*)
      INTEGER, INTENT(IN) :: NC
      REAL(DP) :: DXY, HILF
      INTEGER :: I, IANF, NBEG, J
C
      NBEG=1
      IANF=3
1     CONTINUE
      DO 2 I=IANF+2,NC,2
        DXY=ABS(XY(IANF)-XY(I))+ABS(XY(IANF+1)-XY(I+1))
        IF (DXY.LT.1.E-3) THEN
          IF (I-IANF.GT.2) THEN
            HILF=XY(IANF+2)
            XY(IANF+2)=XY(I)
            XY(I)=HILF
            HILF=XY(IANF+3)
            XY(IANF+3)=XY(I+1)
            XY(I+1)=HILF
            IF (I-IANF.GT.4) THEN
              IF (MOD((I-IANF)/2,2).EQ.1) THEN
                J=I+2
              ELSE
                J=I-2
              ENDIF
              HILF=XY(IANF+4)
              XY(IANF+4)=XY(J)
              XY(J)=HILF
              HILF=XY(IANF+5)
              XY(IANF+5)=XY(J+1)
              XY(J+1)=HILF
            ENDIF
          ENDIF
          IANF=IANF+4
          GOTO 1
        ENDIF
2     CONTINUE
C
      CALL GRJMP (REAL(XY(NBEG),KIND(1.E0)),REAL(XY(NBEG+1),KIND(1.E0)))
      DO 5 I=NBEG+2,IANF,4
5       CALL GRDRW (REAL(XY(I),KIND(1.E0)),REAL(XY(I+1),KIND(1.E0)))
C
      IF (IANF.LT.NC-1) THEN
        NBEG=IANF+2
        IANF=IANF+4
        GOTO 1
      ENDIF
C
      RETURN
      END
