C
C
      SUBROUTINE EIRENE_MUELAM (X1,Y1,VX,VY,X2,Y2,X3,Y3,XL,XM)
 
      USE EIRMOD_PRECISION
 
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1, Y1, VX, VY, X2, Y2, X3, Y3
      REAL(DP), INTENT(OUT) :: XL, XM
      REAL(DP) :: D
 
      D=VX*(Y2-Y3)-VY*(X2-X3)+1.E-20
      XL=((Y2-Y3)*(X2-X1)-(X2-X3)*(Y2-Y1))/D
      XM=(-VY*(X2-X1)+VX*(Y2-Y1))/D
 
      RETURN
      END
