C
C
      SUBROUTINE SETREF(AFF,AFFI,IFLAG,CC1,CC2,CC3)
C
C  SET REFLECTION MATRIX AFF AND INVERSE REFLECTION MATRIX AFFI = AFF
C  INPUT:
C  IFLAG=1:  REFLECTION HYPERPLANE NORMAL VECTOR C1,C2,C3
C            I.E. REFLECTION AT PLANE X*C1+Y*C2+Z*C3+0 (!!!)=0
C
      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: AFF(3,3),AFFI(3,3)
      REAL(DP), INTENT(IN) :: CC1, CC2, CC3
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: C, C1, C2, C3, CN
      INTEGER :: I, J

      C1=CC1
      C2=CC2
      C3=CC3
      IF (IFLAG.EQ.1) THEN
C  NORMALIZE ROTATION AXIS
        C=C1*C1+C2*C2+C3*C3
        IF (C.LE.0.) THEN
          WRITE (6,*) 'WARNING: INVALID REFLEC. PLANE IN SUBR. SETREF'
          WRITE (6,*) 'NO REFLECTION CARRIED OUT'
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
        AFF(1,1)=1.-2.*C1*C1
        AFF(2,1)=  -2.*C2*C1
        AFF(3,1)=  -2.*C3*C1
        AFF(1,2)=  -2.*C1*C2
        AFF(2,2)=1.-2.*C2*C2
        AFF(3,2)=  -2.*C3*C2
        AFF(1,3)=  -2.*C1*C3
        AFF(2,3)=  -2.*C2*C3
        AFF(3,3)=1.-2.*C3*C3
        DO 10 I=1,3
          DO 10 J=1,3
            AFFI(I,J)=AFF(I,J)
10      CONTINUE
      ELSE
      ENDIF
      RETURN
      END
