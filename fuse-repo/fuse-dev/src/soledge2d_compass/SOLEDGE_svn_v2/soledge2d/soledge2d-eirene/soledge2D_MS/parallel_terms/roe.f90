SUBROUTINE ROE(Nx,Nz,LM,LP,FM,FP,OMM,OMP,EIGENM,EIGENP)
  IMPLICIT NONE
  INTEGER*4 :: Nx,Nz
  REAL*8,DIMENSION(1:Nx,0:Nz,1:3),intent(in) :: LM,LP,FM,FP,OMM,OMP
  real*8,dimension(1:Nx,0:Nz,1:3),intent(out):: EIGENM,EIGENP
  INTEGER*4 i,j,k
  REAL*8 ALPHA
  DO k=1,3
     DO I=1,Nx
        do j=0,Nz
           IF(LM(i,j,k)*LP(i,j,k).GT.0.d0) THEN
              ALPHA=ABS(LP(i,j,k))/LP(i,j,k)
              EIGENM(i,j,k) = 0.5d0*(1.d0+ALPHA)*FM(i,j,k)
              EIGENP(i,j,k) = 0.5d0*(1.d0-ALPHA)*FP(i,j,k)
           ELSE
              ALPHA =MAX(ABS(LM(i,j,k)),ABS(LP(i,j,k)))
              EIGENM(i,j,k)=0.5d0*(FM(i,j,k)+ALPHA*OMM(i,j,k))
              EIGENP(i,j,k)=0.5d0*(FP(i,j,k)-ALPHA*OMP(i,j,k))
           ENDIF
        end do
     ENDDO
  ENDDO
END SUBROUTINE ROE
