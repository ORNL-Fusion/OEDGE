SUBROUTINE ROE_e(Nx,Nz,LM,LP,FM,FP,OMM,OMP,EIGENM,EIGENP)
  IMPLICIT NONE
  INTEGER*4 :: Nx,Nz
  REAL*8,DIMENSION(1:Nx,0:Nz),intent(in) :: LM,LP,FM,FP,OMM,OMP
  real*8,dimension(1:Nx,0:Nz),intent(out):: EIGENM,EIGENP
  INTEGER*4 i,j,k
  REAL*8 ALPHA
  DO I=1,Nx
     do j=0,Nz
        IF(LM(i,j)*LP(i,j).GT.0.d0) THEN
           ALPHA=ABS(LP(i,j))/LP(i,j)
           EIGENM(i,j) = 0.5d0*(1.d0+ALPHA)*FM(i,j)
           EIGENP(i,j) = 0.5d0*(1.d0-ALPHA)*FP(i,j)
        ELSE
           ALPHA =MAX(ABS(LM(i,j)),ABS(LP(i,j)))
           EIGENM(i,j)=0.5d0*(FM(i,j)+ALPHA*OMM(i,j))
           EIGENP(i,j)=0.5d0*(FP(i,j)-ALPHA*OMP(i,j))
        ENDIF
     end do
  ENDDO
END SUBROUTINE ROE_E
