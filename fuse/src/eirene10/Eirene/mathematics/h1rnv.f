*
      SUBROUTINE H1RNV(RVEC,LEN)
*
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: RVEC(1)
      INTEGER, INTENT(IN) :: LEN
      REAL(DP) :: UNI
      INTEGER :: IVEC
      CHARACTER(16) :: FLAG,CHECK
      REAL(DP) :: U, C, CD, CM
      INTEGER :: I, J
      COMMON /RASET1/ U(97),C,CD,CM,I,J
      COMMON /RASET2/ FLAG
*
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      DATA CHECK /'H1RN INITIALISED'/
*
      IF (FIRST) THEN
         IF (FLAG .NE. CHECK) THEN
            WRITE(iunout,*) 
     >                 ' H1RNV (RANMAR): CALL H1RNIN OR H1RNIV BEFORE',
     >                 ' CALLING H1RN.'
            STOP
         ELSE
            FIRST = .FALSE.
         ENDIF
      ENDIF
*
      DO 200 IVEC = 1,LEN
  190    CONTINUE
         UNI = U(I)-U(J)
         IF(UNI .LT. 0.) UNI = UNI + 1.
         U(I) = UNI
         I = I - 1
         IF( I .EQ. 0) I=97
         J = J - 1
         IF( J .EQ. 0) J=97
         C = C - CD
         IF( C .LT. 0) C = C + CM
         UNI = UNI-C
         IF(UNI .LE. 0.) UNI = UNI + 1.
         IF(UNI .GE. 1.) GOTO 190
         RVEC(IVEC) = UNI
  200 CONTINUE
      RETURN

C     The following ENTRY is for reinitialization of EIRENE
      
      ENTRY H1RNV_REINIT
      FIRST = .TRUE.
      FLAG = 'H1RN UNINITIALIS'
      return
      END
