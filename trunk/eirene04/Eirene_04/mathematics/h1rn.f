C
*
      FUNCTION H1RN(DUMMY)
*
*#**********************************************************************
*# RANDOM NUMBER GENERATOR AS ADVOCATED BY F. JAMES FROM PROPOSAL OF   *
*# MARSAGLIA AND ZAMAN FSU-SCRI-87-50 AND MODIFIED BY F. JAMES 1988 TO *
*# PRODUCE VECTOR OF NUMBERS.                                          *
*# ENTRIES ARE:                                                        *
*#     FUNCTION    H1RN(DUMMY)     SINGLE RANDOM NUMBER                *
*#     SUBROUTINE  H1RNV(VEC,LEN)  VECTOR OF RANDOM NUMBERS            *
*#     SUBROUTINE  H1RNIN(IJ,KL)   INITIALISE WITH SEEDS               *
*#     SUBROUTINE  H1RNIV(VEC)     INITIALISE/RESTART WITH SEED ARRAY  *
*#     SUBROUTINE  H1RNSV(VEC)     SAVE SEED ARRAY VEC(100)            *
*#                                                                     *
*# NOTE: -H1RNIN OR H1RNIV MUST BE CALLED BEFORE GENERATING ANY        *
*#        RANDOM NUMBER(S).                                            *
*#       -H1RNSV SAVES SEED ARRAY INTO VEC(100) ONLY. THE USER HAS TO  *
*#        OUTPUT IT.                                                   *
*#                                                                     *
*# CHANGED BY: G. GRINDHAMMER AT: 90/03/14                             *
*# REASON :                                                            *
*#**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: DUMMY
      REAL(DP) :: H1RN
      INTEGER :: ISEED1, ISEED2
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
cpb         WRITE(6,*) ' H1RN (RANMAR): INITIALIZED WITH DEFAULT SEED'
            ISEED1      = 12345
            ISEED2      = 98765
            CALL H1RNIN(ISEED1,ISEED2)
            FIRST = .FALSE.
         ELSE
            FIRST = .FALSE.
         ENDIF
      ENDIF
*
  100 CONTINUE
      H1RN = U(I)-U(J)
      IF(H1RN .LT. 0.) H1RN = H1RN + 1.
      U(I) = H1RN
      I = I - 1
      IF( I .EQ. 0) I=97
      J = J - 1
      IF( J .EQ. 0) J=97
      C = C - CD
      IF( C .LT. 0) C = C + CM
      H1RN = H1RN-C
      IF(H1RN .LE. 0.) H1RN = H1RN + 1.
      IF(H1RN .GE. 1.) GOTO 100
      RETURN
      END
