CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ****s* INTERPOLATION/BILINEAR_INT
C NAME
C     BILINEAR_INT (D, RX, RY, Z)
C DESCRIPTION
C     Bilinear interpolation of data
C INPUTS
C     REAL*8, DIMENSION(2,2):: D  ! data to interpolate
C     REAL*8    :: RX, RY ! relative position at which interpolations should be done
C OUTPUT
C     REAL*8    :: Z  ! interpolated data
C     ******
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE EIRENE_BILINEAR_INT (D, RX, RY, Z)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
 
!---- input variables
      REAL(DP), DIMENSION(2,2), INTENT(IN):: D
      REAL(DP), INTENT(IN)  :: RX, RY
!---- output variables
      REAL(DP), INTENT(OUT)  :: Z
 
!---- local variables
      REAL(DP)  :: Z1, Z2
 
      Z1 = D(1,1)*RX + D(2,1)*(1-RX)
      Z2 = D(1,2)*RX + D(2,2)*(1-RX)
      Z  = Z1*RY + Z2*(1-RY)
 
      RETURN
      END SUBROUTINE EIRENE_BILINEAR_INT
