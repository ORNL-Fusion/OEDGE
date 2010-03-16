

      SUBROUTINE INSERT_POINT (BAUM,X,Y,Z,DIST,IC)
C
      USE PRECISION
      USE PARMMOD
      USE CTETRA
      USE COMPRT, ONLY: IUNOUT
      USE module_avltree

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X,Y,Z,DIST
      INTEGER, INTENT(OUT) :: IC
      INTEGER I, ICOOR
      type(TAVLTree), pointer :: baum
      logical :: inserted

      IC = 0
      IC=NCOOR+1
      inserted=.false.
      call insert (baum, x, y, z, dist, ic, inserted)
      IF (.not.inserted) RETURN

      NCOOR=NCOOR+1
      IF (NCOOR > NCOORD) THEN
        WRITE (iunout,*) ' ALLOWED NUMBER OF COORDINATES EXCEEDED '
        WRITE (iunout,*) ' INCREASE NCOORD '
        CALL EXIT_OWN(1)
      END IF
      XTETRA(NCOOR) = X
      YTETRA(NCOOR) = Y
      ZTETRA(NCOOR) = Z
      IC = NCOOR
      RETURN
      END



