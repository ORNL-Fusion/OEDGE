C
C
C*DK MARPRM
      SUBROUTINE MASPRM (A,NA,MA,B,NB,MB,IERR)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A, B
      INTEGER, INTENT(IN) :: NA, MA, NB, MB
      INTEGER, INTENT(INOUT) :: IERR
      IERR=IERR+1
      WRITE (iunout,*) 'PARAMETER ERROR DETECTED '
      WRITE (iunout,*) A(1:NA),' MUST BE >= ',B(1:NB)
      WRITE (iunout,*) A(1:NA),' = ',MA
      WRITE (iunout,*) B(1:NB),' = ',MB
      RETURN
      END
