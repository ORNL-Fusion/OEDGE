C
C
      SUBROUTINE EIRENE_SETROT(AFF,AFFI,IFLAG,CC1,CC2,CC3,CC4)
C
C  SET ROTATION MATRIX AFF AND INVERSE ROTATION MATRIX AFFI
C  INPUT:
C  IFLAG=1:  ROTATION AXIS C1,C2,C3 AND ROTATION ANGLE C4 DEGREES
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
 
      REAL(DP), INTENT(OUT) :: AFF(3,3),AFFI(3,3)
      REAL(DP), INTENT(IN) :: CC1, CC2, CC3, CC4
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: C, C1, C2, C3, C4, CAL, SAL, CN, ANG, PI
      INTEGER :: I, J
 
      DATA PI/3.141592654/
 
      C1=CC1
      C2=CC2
      C3=CC3
      C4=CC4
      IF (IFLAG.EQ.1) THEN
C  NORMALIZE ROTATION AXIS
        C=C1*C1+C2*C2+C3*C3
        IF (C.LE.0.) THEN
          WRITE (iunout,*)
     .    'WARNING: INVALID ROTATION AXIS IN SUBR. SETROT'
          WRITE (iunout,*) 'NO ROTATION CARRIED OUT'
          DO 1 J=1,3
            DO 1 I=1,3
              AFF(I,J)=0.
              AFFI(I,J)=0.
1         CONTINUE
          DO 2 J=1,3
            AFF(J,J)=1.
            AFFI(J,J)=1.
2         CONTINUE
          RETURN
        ENDIF
        CN=SQRT(C)
        C1=C1/CN
        C2=C2/CN
        C3=C3/CN
        ANG=C4*PI/180.
        CAL=COS(ANG)
        SAL=SIN(ANG)
        AFF(1,1)=CAL+(1-CAL)*C1*C1
        AFF(2,1)=    (1-CAL)*C2*C1+SAL*C3
        AFF(3,1)=    (1-CAL)*C3*C1-SAL*C2
        AFF(1,2)=    (1-CAL)*C1*C2-SAL*C3
        AFF(2,2)=CAL+(1-CAL)*C2*C2
        AFF(3,2)=    (1-CAL)*C3*C2+SAL*C1
        AFF(1,3)=    (1-CAL)*C1*C3+SAL*C2
        AFF(2,3)=    (1-CAL)*C2*C3-SAL*C1
        AFF(3,3)=CAL+(1-CAL)*C3*C3
        CAL=COS(-ANG)
        SAL=SIN(-ANG)
        AFFI(1,1)=CAL+(1-CAL)*C1*C1
        AFFI(2,1)=    (1-CAL)*C2*C1+SAL*C3
        AFFI(3,1)=    (1-CAL)*C3*C1-SAL*C2
        AFFI(1,2)=    (1-CAL)*C1*C2-SAL*C3
        AFFI(2,2)=CAL+(1-CAL)*C2*C2
        AFFI(3,2)=    (1-CAL)*C3*C2+SAL*C1
        AFFI(1,3)=    (1-CAL)*C1*C3+SAL*C2
        AFFI(2,3)=    (1-CAL)*C2*C3-SAL*C1
        AFFI(3,3)=CAL+(1-CAL)*C3*C3
      ELSE
      ENDIF
      RETURN
      END
