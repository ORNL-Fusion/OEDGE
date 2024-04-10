C
C
      SUBROUTINE EIRENE_PSIDE (G,R,P11,P12,P21,P22,L1,L2,EPS)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CRECH
      USE EIRMOD_CPLOT
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: G(4,2), R(4,2)
      REAL(DP), INTENT(IN) :: P11, P12, P21, P22, EPS
      LOGICAL, INTENT(IN) :: L1, L2
      REAL(DP) :: XLA, XP, YP, XMU
      INTEGER :: I, IS
 
      IS=0
      DO 100 I=1,4
      CALL EIRENE_GSP
     .  (G(I,1),G(I,2),R(I,1),R(I,2),P11,P12,P21-P11,P22-P12,
     .          XLA,XMU,EPS)
      IF (XLA.GE.0..AND.XLA.LE.1..AND.XMU.GE.0..AND.XMU.LE.1.) THEN
        XP=G(I,1)+XMU*R(I,1)
        YP=G(I,2)+XMU*R(I,2)
        IF (L1) THEN
          IF (LZR) THEN
            CALL GRJMP
     .  (REAL(P11,KIND(1.E0)),REAL(P12,KIND(1.E0)))
            CALL GRDRW (REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
          ENDIF
          IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
            CALL EIRENE_STCOOR (P11,P12,0)
            CALL EIRENE_STCOOR (XP,YP,1)
          ENDIF
          RETURN
        ELSE IF (L2) THEN
          IF (LZR) THEN
            CALL GRJMP (REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
            CALL GRDRW
     .  (REAL(P21,KIND(1.E0)),REAL(P22,KIND(1.E0)))
          ENDIF
          IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
            CALL EIRENE_STCOOR (XP,YP,0)
            CALL EIRENE_STCOOR (P21,P22,1)
          ENDIF
          RETURN
        ELSE
          IS=IS+1
          IF (IS.EQ.1) THEN
            IF (LZR) CALL GRJMP (REAL(XP,KIND(1.E0)),
     .                           REAL(YP,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL EIRENE_STCOOR
     .  (XP,YP,0)
          ELSEIF (IS.EQ.2) THEN
            IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                           REAL(YP,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL EIRENE_STCOOR
     .  (XP,YP,1)
            RETURN
          ENDIF
        ENDIF
      ENDIF
100   CONTINUE
      RETURN
      END
