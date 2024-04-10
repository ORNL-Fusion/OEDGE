subroutine set_test_geometry()
  use test_var
  use all_variables, only : global_parameters, zones
  use Mphysics
  implicit none
  integer*4 :: i,j,Nx,Nz,k
  real*8 :: r,theta
  real*8 :: test_theta,test_r
  call allocate_geometry()
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           r=test_r(zones(k)%mesh%x(i,j))
           theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))
           zones(k)%mesh%Rgeom(i,j)=test_R0+r*cos(theta)
           zones(k)%mesh%Zgeom(i,j)=r*sin(theta)
           zones(k)%mesh%Br(i,j)=test_psi0/(test_a*zones(k)%mesh%Rgeom(i,j))*sin(theta)
           zones(k)%mesh%Bz(i,j)=-test_psi0/(test_a*zones(k)%mesh%Rgeom(i,j))*cos(theta)
           zones(k)%mesh%Bphi(i,j)=test_B0*test_R0/zones(k)%mesh%Rgeom(i,j)
           r=test_r(zones(k)%mesh%x_minus_1half(i,j))
           theta=test_theta(zones(k)%mesh%z_minus_1half(i,j),zones(k)%mesh%x_minus_1half(i,j))
           zones(k)%mesh%Rcorner(i,j)=test_R0+r*cos(theta)
           zones(k)%mesh%Zcorner(i,j)=r*sin(theta)
        end do
        !Nz+1
        r=test_r(zones(k)%mesh%x_minus_1half(i,Nz))
        theta=test_theta(zones(k)%mesh%z_plus_1half(i,Nz),zones(k)%mesh%x_minus_1half(i,Nz))
        zones(k)%mesh%Rcorner(i,Nz+1)=test_R0+r*cos(theta)
        zones(k)%mesh%Zcorner(i,Nz+1)=r*sin(theta)
     end do
     !Nx+1
     do j=1,Nz
        r=test_r(zones(k)%mesh%x_plus_1half(Nx,j))
        theta=test_theta(zones(k)%mesh%z_minus_1half(Nx,j),zones(k)%mesh%x_plus_1half(Nx,j))
        zones(k)%mesh%Rcorner(Nx+1,j)=test_R0+r*cos(theta)
        zones(k)%mesh%Zcorner(Nx+1,j)=r*sin(theta)
     end do
     !Nx+1 Nz+1
     r=test_r(zones(k)%mesh%x_plus_1half(Nx,Nz))
     theta=test_theta(zones(k)%mesh%z_plus_1half(Nx,Nz),zones(k)%mesh%x_plus_1half(Nx,Nz))
     zones(k)%mesh%Rcorner(Nx+1,Nz+1)=test_R0+r*cos(theta)
     zones(k)%mesh%Zcorner(Nx+1,Nz+1)=r*sin(theta)
  end do
  !Neighboring points where needed
  do j=1,Nz
     !zone 1 south
     r=test_r(zones(1)%mesh%x_minus_1half(1,j))
     theta=test_theta(zones(1)%mesh%z(1,j),zones(1)%mesh%x_minus_1half(1,j))
     zones(1)%mesh%Rgeom(0,j)=test_R0+r*cos(theta)
     zones(1)%mesh%Zgeom(0,j)=r*sin(theta)
     zones(1)%mesh%Br(0,j)=test_psi0/(test_a*zones(1)%mesh%Rgeom(0,j))*sin(theta)
     zones(1)%mesh%Bz(0,j)=-test_psi0/(test_a*zones(1)%mesh%Rgeom(0,j))*cos(theta)
     zones(1)%mesh%Bphi(0,j)=test_B0*test_R0/zones(1)%mesh%Rgeom(0,j)
     !zone 2 north
     r=test_r(zones(2)%mesh%x_plus_1half(Nx,j))
     theta=test_theta(zones(2)%mesh%z(Nx,j),zones(2)%mesh%x_plus_1half(Nx,j))
     zones(2)%mesh%Rgeom(Nx+1,j)=test_R0+r*cos(theta)
     zones(2)%mesh%Zgeom(Nx+1,j)=r*sin(theta)
     zones(2)%mesh%Br(Nx+1,j)=test_psi0/(test_a*zones(2)%mesh%Rgeom(Nx+1,j))*sin(theta)
     zones(2)%mesh%Bz(Nx+1,j)=-test_psi0/(test_a*zones(2)%mesh%Rgeom(Nx+1,j))*cos(theta)
     zones(2)%mesh%Bphi(Nx+1,j)=test_B0*test_R0/zones(2)%mesh%Rgeom(Nx+1,j)
  end do
  do k=1,2
     zones(k)%mesh%B=sqrt(Zones(k)%mesh%Br**2.D0+Zones(k)%mesh%Bz**2.D0+Zones(k)%mesh%Bphi**2.D0)
  end do
  call save_mesh_and_geometry()
end subroutine set_test_geometry
