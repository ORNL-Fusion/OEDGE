subroutine compute_weno_coef()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 i,j,k,Nx,Nz
  real*8,allocatable:: zp12(:,:),zm12(:,:)
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     allocate(Zones(k)%Weno%C20_i_xp12(1:Nx,0:Nz))
     allocate(Zones(k)%Weno%C21_i_xp12(1:Nx,0:Nz))
     allocate(Zones(k)%Weno%C20_ip1_xp12(1:Nx,0:Nz))
     allocate(Zones(k)%Weno%C21_ip1_xp12(1:Nx,0:Nz))
     allocate(zp12(1:Nx,-1:Nz+2))
     zp12(1:Nx,-1:Nz+2)=Zones(k)%mesh%z_plus_1half(1:Nx,-1:Nz+2)
     allocate(zm12(1:Nx,-1:Nz+2))
     zm12(1:Nx,-1:Nz+2)=Zones(k)%mesh%z_minus_1half(1:Nx,-1:Nz+2)
     do i=1,Nx
        do j=0,Nz
           !polynome sur Ii en x_i+1/2
           Zones(k)%Weno%C20_i_xp12(i,j)=(zp12(i,j+1)-zp12(i,j))/(zp12(i,j+1)-zm12(i,j-1))
           Zones(k)%Weno%C21_i_xp12(i,j)=(zm12(i,j-1)-zp12(i,j))/(zm12(i,j-1)-zp12(i,j+1))
           !polynome sur Ii+1 en x_i+1/2
           Zones(k)%Weno%C20_ip1_xp12(i,j)=(zp12(i,j+2)-zp12(i,j))/(zp12(i,j+2)-zm12(i,j))
           Zones(k)%Weno%C21_ip1_xp12(i,j)=(zp12(i,j)-zm12(i,j))/(zp12(i,j+2)-zm12(i,j))
        end do
     end do
     deallocate(zp12,zm12)
  end do
end subroutine compute_weno_coef

subroutine weno2fv_nonugrid(Zone,Nx,Nz,var,u_xp12_m,u_xp12_p)
  use MZone
  implicit none
  integer*4 Nx,Nz
  Type(TZone) :: Zone
  REAL*8,DIMENSION(1:Nx,-1:Nz+2) :: var
  REAL*8,DIMENSION(1:Nx,0:Nz) :: u_xp12_m,u_xp12_p
  integer*4 i,j
  real*8 :: alpha20,IS20,omega20
  real*8 :: alpha21,IS21,omega21
  real*8 :: b00,b01,b10,b11
  real*8,dimension(1:Nx,-1:Nz+2) :: zp12,zm12
  real*8 epsilon
  epsilon=1.D-6
  zp12(1:Nx,-1:Nz+2)=Zone%mesh%z_plus_1half(1:Nx,-1:Nz+2)
  zm12(1:Nx,-1:Nz+2)=Zone%mesh%z_minus_1half(1:Nx,-1:Nz+2)
  do i=1,Nx
     do j=0,Nz
        !u_xp12_m (polynom p_i)
        b00=var(i,j-1)+(2.d0*zp12(i,j)-zm12(i,j-1)-zm12(i,j))*&
             (var(i,j)-var(i,j-1))/(zp12(i,j)-zm12(i,j-1))
        b01=2*(var(i,j)-var(i,j-1))/(zp12(i,j)-zm12(i,j-1))
        b11=2*(var(i,j+1)-var(i,j))/(zp12(i,j+1)-zm12(i,j))
        b10=var(i,j)+(2.d0*zp12(i,j)-zm12(i,j)-zm12(i,j+1))*&
             (var(i,j+1)-var(i,j))/(zp12(i,j+1)-zm12(i,j))
        IS20=(b01*(zp12(i,j)-zm12(i,j)))**2.d0
        IS21=(b11*(zp12(i,j)-zm12(i,j)))**2.d0
        alpha20=(zone%Weno%C20_i_xp12(i,j)/(epsilon+IS20)**2.d0)*2./3.
        alpha21=(zone%Weno%C21_i_xp12(i,j)/(epsilon+IS21)**2.d0)*1./3.
        omega20=alpha20/(alpha20+alpha21)
        omega21=alpha21/(alpha20+alpha21)
        u_xp12_m(i,j)=omega20*(b00)+omega21*(b10)
        !u_xp12_p (polynome p_i+1)
        b00=var(i,j)+(2.d0*zp12(i,j+1)-zm12(i,j)-zm12(i,j+1))*&
             (var(i,j+1)-var(i,j))/(zp12(i,j+1)-zm12(i,j))
        b01=2*(var(i,j+1)-var(i,j))/(zp12(i,j+1)-zm12(i,j))
        b11=2*(var(i,j+2)-var(i,j+1))/(zp12(i,j+2)-zm12(i,j+1))
        b10=var(i,j+1)+(2.d0*zp12(i,j+1)-zm12(i,j+1)-zm12(i,j+2))*&
             (var(i,j+2)-var(i,j+1))/(zp12(i,j+2)-zm12(i,j+1))
        IS20=(b01*(zp12(i,j+1)-zm12(i,j+1)))**2.d0
        IS21=(b11*(zp12(i,j+1)-zm12(i,j+1)))**2.d0
        alpha20=(zone%Weno%C20_ip1_xp12(i,j)/(epsilon+IS20)**2.d0)*1./3.
        alpha21=(zone%Weno%C21_ip1_xp12(i,j)/(epsilon+IS21)**2.d0)*2./3.
        omega20=alpha20/(alpha20+alpha21)
        omega21=alpha21/(alpha20+alpha21)
        u_xp12_p(i,j)=omega20*(b00+b01*(zp12(i,j)-zp12(i,j+1)))&
             +omega21*(b10+b11*(zp12(i,j)-zp12(i,j+1)))
     end do
  end do
end subroutine weno2fv_nonugrid

subroutine weno2_unif(Nx,Nz,var,u_xp12_m,u_xp12_p)
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
  epsilon=1.D-7
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
