subroutine compute_weno_coef()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: i,j,k,Nx,Nz
  integer*4 :: r,j_,m,l,q
  real*8 :: t1,t2,t3
  real*8,allocatable :: zp12(:,:),zm12(:,:)
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     allocate(Zones(k)%Weno(1:Nx,0:Nz))
     allocate(zp12(1:Nx,-1:Nz+2))
     zp12(1:Nx,-1:Nz+2)=Zones(k)%mesh%z_plus_1half(1:Nx,-1:Nz+2)
     allocate(zm12(1:Nx,-1:Nz+2))
     zm12(1:Nx,-1:Nz+2)=Zones(k)%mesh%z_minus_1half(1:Nx,-1:Nz+2)
     do i=1,Nx
        do j=0,Nz
           zones(k)%Weno(i,j)%C=0.D0
           do r=0,1
              do j_=0,1
                 t1=0.D0
                 do m=j_+1,2
                    t2=0.D0
                    do l=0,2
                       if(l.ne.m) then
                          t3=1.D0
                          do q=0,2
                             if((q.ne.m).and.(q.ne.l)) then
                                t3=t3*(zp12(i,j)-zm12(i,j-r+q))
                             end if
                          end do
                          t2=t2+t3
                       end if
                    end do
                    t3=1.D0
                    do l=0,2
                       if(l.ne.m) then
                          t3=t3*(zm12(i,j-r+m)-zm12(i,j-r+l))
                       end if
                    end do
                    t1=t1+t2/t3
                 end do
                 zones(k)%weno%C(r,j_)=t1*(zp12(i,j-r+j_)-zm12(i,j-r+j_))
              end do
           end do
        end do
     end do

     !Ctilde
     do i=1,Nx
        do j=0,Nz
           zones(k)%Weno(i,j)%Ctilde=0.D0
           do r=0,1
              do j_=0,1
                 t1=0.D0
                 do m=j_+1,2
                    t2=0.D0
                    do l=0,2
                       if(l.ne.m) then
                          t3=1.D0
                          do q=0,2
                             if((q.ne.m).and.(q.ne.l)) then
                                t3=t3*(zp12(i,j)-zp12(i,j-r+q))
                             end if
                          end do
                          t2=t2+t3
                       end if
                    end do
                    t3=1.D0
                    do l=0,2
                       if(l.ne.m) then
                          t3=t3*(zp12(i,j-r+m)-zp12(i,j-r+l))
                       end if
                    end do
                    t1=t1+t2/t3
                 end do
                 zones(k)%weno%Ctilde(r,j_)=t1*(zp12(i,j+1-r+j_)-zm12(i,j+1-r+j_))
              end do
           end do
        end do
     end do

     deallocate(zm12,zp12)
  end do
end subroutine compute_weno_coef


subroutine compute_weno_nonuni(Zone,Nx,Nz,var,u_xp12_m,u_xp12_p)
  use MZone
  implicit none
  integer*4 Nx,Nz
  Type(TZone) :: Zone
  REAL*8,DIMENSION(1:Nx,-1:Nz+2) :: var
  REAL*8,DIMENSION(1:Nx,0:Nz) :: u_xp12_m,u_xp12_p
  integer*4 :: i,j
  integer*4 :: r,j_
  real*8 :: v(0:1)
  real*8 :: beta0,beta1
  real*8 :: alpha0,alpha1
  real*8 :: omega0,omega1
  real*8,parameter :: epsilon=1.d-6
  do i=1,Nx
     do j=0,Nz
        ! var i+1/2 minus
        do r=0,1
           v(r)=0.D0
           do j_=0,1
              v(r)=v(r)+zone%Weno(i,j)%C(r,j_)*var(i,j-r+j_)
           end do
        end do
        beta0=(var(i,j+1)-var(i,j))**2
        beta1=(var(i,j)-var(i,j-1))**2
        alpha0=2./3./(epsilon+beta0)**2
        alpha1=1./3./(epsilon+beta1)**2
        omega0=alpha0/(alpha0+alpha1)
        omega1=alpha1/(alpha0+alpha1)
        u_xp12_m(i,j)=omega0*v(0)+omega1*v(1)
        ! var i+1/2 plus
        do r=0,1
           v(r)=0.D0
           do j_=0,1
              v(r)=v(r)+zone%Weno(i,j)%Ctilde(r,j_)*var(i,j+1-r+j_)
           end do
        end do
        beta0=(var(i,j+2)-var(i,j+1))**2
        beta1=(var(i,j+1)-var(i,j))**2
        alpha0=1./3./(epsilon+beta0)**2
        alpha1=2./3./(epsilon+beta1)**2
        omega0=alpha0/(alpha0+alpha1)
        omega1=alpha1/(alpha0+alpha1)
        u_xp12_p(i,j)=omega0*v(0)+omega1*v(1)
     end do
  end do
end subroutine compute_weno_nonuni
