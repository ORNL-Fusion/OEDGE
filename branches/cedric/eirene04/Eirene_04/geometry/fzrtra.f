C
C
      SUBROUTINE FZRTRA(X,Z,PH,NNEW)

      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGRID

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X, Z
      REAL(DP), INTENT(INOUT) :: PH
      INTEGER, INTENT(OUT) :: NNEW
      REAL(DP) :: Z1, X01
      INTEGER :: LEARCA
C
C
C  FIND X,Z, NNEW,   FROM X,PH   (X=X??)
      Z1=ZSURF(1)
      PH=MOD(PH+PI2A-Z1,PI2A)+Z1
      NNEW=LEARCA(PH,ZSURF,1,NTTRA,1,'FZRTRA ')
      IF (NNEW.LE.0.OR.NNEW.GT.NTTRAM) THEN
        WRITE (6,*) 'NT OUT OF RANGE IN FZRTRA '
        WRITE (6,*) PH,ZHALF,NNEW
        CALL EXIT_OWN(1)
      ENDIF
      X01=X+RMTOR
      CALL FZRTRI(X,Z,NNEW,X01,PH,NNEW)
      RETURN
      END