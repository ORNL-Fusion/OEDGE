subroutine weno2_unif(Nx,Nz,var,u_xp12_m,u_xp12_p)
  use all_variables, only : global_variables
  implicit none
  integer*4 Nx,Nz
  REAL*8,DIMENSION(1:Nx,-1:Nz+2) :: var
  REAL*8,DIMENSION(1:Nx,0:Nz) :: u_xp12_m,u_xp12_p
  integer*4 i,j,k
  real*8  utemp(4),dif1(4)
  INTEGEr*4 i0,i1,i2
  real*8 px,ind1,ind2
  real*8 fac1,fac2
  real*8 epsilon
  epsilon=global_variables%epsilon_weno
  do i=1,Nx
     do j=0,Nz
        !u_xp12_m (polynom p_i)
        utemp(1) = var(i,j-1)
        utemp(2) = var(i,j)
        utemp(3) = var(i,j+1)
        utemp(4) = var(i,j+2)

        DO k=1,3
           dif1(k) = utemp(k+1) - utemp(k)
        ENDDO

        i0=3
        i0=i0-1
        i1=i0+1
        i2=i0+2

        ind1=(dif1(2)**2)
        ind2=(dif1(3)**2)

        ind1=1.d0/(epsilon+ind1)**2
        ind2=1.d0/(epsilon+ind2)**2/2.d0
        fac1=ind1/(ind1+ind2)
        fac2=ind2/(ind1+ind2)

        u_xp12_p(i,j)=   utemp(i0) * ( 1.d0*fac1          )*0.5d0        &
             +utemp(i1) * ( 1.d0*fac1+3.d0*fac2)*0.5d0        &
             +utemp(i2) * (-1.d0*fac2          )*0.5d0

        i0=2
        i0=i0-1
        i1=i0+1
        i2=i0+2
        ind1=(dif1(1)**2)
        ind2=(dif1(2)**2)

        ind1=1.d0/(epsilon+ind1)**2/2.d0
        ind2=1.d0/(epsilon+ind2)**2

        fac1=ind1/(ind1+ind2)
        fac2=ind2/(ind1+ind2)

        u_xp12_m(i,j)=   utemp(i0) * (-1.d0*fac1          )*0.5d0        &
             +utemp(i1) * ( 3.d0*fac1+1.d0*fac2)*0.5d0        &
             +utemp(i2) * ( 1.d0*fac2          )*0.5d0
     end do
  end do
end subroutine weno2_unif
