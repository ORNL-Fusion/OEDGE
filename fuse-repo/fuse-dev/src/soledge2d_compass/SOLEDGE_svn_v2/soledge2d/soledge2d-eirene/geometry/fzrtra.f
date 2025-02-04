!pb  05.10.06 : intent of arguments X and Z corrected
C
C
      SUBROUTINE EIRENE_FZRTRA(X,Z,PH,NNEW)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(INOUT) :: X, Z
      REAL(DP), INTENT(INOUT) :: PH
      INTEGER, INTENT(OUT) :: NNEW
      REAL(DP) :: Z1, X01,XX
      INTEGER :: EIRENE_LEARCA
C
C
C  FIND X,Z, NNEW,   FROM X,PH   (X=X??)
      Z1=ZSURF(1)
      PH=MOD(PH+PI2A-Z1,PI2A)+Z1
      NNEW=EIRENE_LEARCA(PH,ZSURF,1,NTTRA,1,'FZRTRA ')
      IF (NNEW.LE.0.OR.NNEW.GT.NTTRAM) THEN
        WRITE (iunout,*) 'NT OUT OF RANGE IN FZRTRA '
        WRITE (iunout,*) PH,ZHALF,NNEW
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      X01=X+RMTOR
      CALL EIRENE_FZRTRI(X,Z,NNEW,X01,PH,NNEW)
      RETURN
      END
