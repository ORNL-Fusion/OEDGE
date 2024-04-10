SUBROUTINE PROJECTION_p(Nx,Nz,VEC,VARIN,VAROUT)
  IMPLICIT NONE
  INTEGER*4 :: Nx,Nz
  REAL*8,DIMENSION(0:Nx,1:Nz,1:3,1:3) :: VEC
  REAL*8,DIMENSION(0:Nx,1:Nz,1:3) :: VARIN,VAROUT
  INTEGER*4 :: i,j,k,n
  VAROUT(:,:,:) = 0.d0
  DO k=1,3
     DO n=1,3
        DO I=0,Nx
           do j=1,Nz
              VAROUT(i,j,n)=VAROUT(i,j,n)+VEC(i,j,n,k)*VARIN(i,j,k)
           end do
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE PROJECTION_P
